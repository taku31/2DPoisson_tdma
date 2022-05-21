% 初期化
clear all;

% パラメーター
Lx = 1;% 計算領域
Ly = 1;% 計算領域
nx = 20;% 要素数
ny = 20;% 要素数
dLx = Lx / nx;% 要素サイズ
dLy = Ly / ny;% 要素サイズ
k = 1;% 熱伝導率

% 座標の生成
x = 0 : Lx / nx : Lx;
y = Ly : - Ly / ny : 0;
[X, Y] = meshgrid(x, y);

% 配列確保
T2d = zeros(nx, ny);

% 熱源の設定
q2d = zeros(nx, ny);
%q2d(round(nx / 2), round(ny / 2)) = 1000;% 中央に単位熱源を設置する

% 時間計測開始
tic

% 連立方程式の解法を選択
method = 'TDMA';

switch method

    case 'SOR'

        T2d = SOR(T2d, q2d, k, nx, ny, dLx, dLy);

    case 'TDMA'

        alpha = 1 / (dLx)^2;
        beta = 1 / (dLy)^2;
        gamma = - 2 * (alpha + beta);
        a = zeros(1, nx - 1);
        b = zeros(nx * ny, 1);
        c = zeros(1, nx - 1);
        d = zeros(1, nx);
        T2d_new = zeros(nx, ny);

        % 反復法のパラメーター
        eps = 10^(-14);% 許容誤差
        ite_max = nx * ny *20;% 最大反復回数
        omega = 1.3;% 緩和係数

        for ite = 1 : ite_max

            error = 0;

            %ite

            for i = 2 : nx - 1

                % ±x側境界条件
                T2d(1, :) = 1;% ディリクレ条件
                T2d(end, :) = T2d(end - 1, :) ; % ノイマン条件
                %T2d(end, :) = 0; % ディリクレ条件

                % j = 1 の時 (-y側の境界条件)
                b(1) = 1;
                c(1) = 0;% ディリクレ条件
                % c(1) = -1;% ノイマン条件
                d(1) = 0;

                % j = 2 ~ ny - 1 の時(内部領域)
                for j = 2 : ny - 1
                    a(j - 1) = beta;
                    b(j) = gamma / omega;
                    c(j) = beta;
                    d(j) = - alpha * T2d(i + 1, j) - alpha * T2d(i - 1, j) - q2d(i, j) / k + (1 - omega) * gamma * T2d(i, j)/ omega;
                end

                % j = ny の時(+y側の境界条件)
                a(ny - 1) = 0;% ディリクレ条件
                % a(ny - 1) = -1;% ノイマン条件
                b(ny) = 1;
                d(ny) = 0;

                % TDMA
                T2d_new(i, :) = TDMA(a, b, c, d);

                % 誤差の計算
                for j = 2 : ny - 1
                    cor = T2d(i, j) - T2d_new(i, j);
                    error = max(error, abs(cor));
                    T2d(i, j) = T2d_new(i, j);
                end

                % 角部分は境界条件の平均とする
                T2d(1, 1) = (T2d(1, 2) + T2d(2, 1)) / 2;
                T2d(1, end) = (T2d(1, end - 1) + T2d(2, end)) / 2;
                T2d(end, 1) = (T2d(end, 2) + T2d(end - 1, 1)) / 2;
                T2d(end, end) = (T2d(end - 1, end) + T2d(end, end - 1)) / 2;

            end

            if error < eps %&& ite> 100 % 収束条件が満たされたらループを抜ける。
                break
            end

            if ite >= ite_max
                disp('最大反復回数に達しました。収束条件を満たしていません。');
            end

        end

    otherwise

        error("選択された手法：" + method + "はサポートしていないです。" + ...
            "SOR,TDMAから選択してください。");

end

% 時間計測終了
toc

% 可視化
vis_contour(T2d, x, y, Lx, Ly, 2)

%% 以下関数
function[T2d] = SOR(T2d, q2d, k, nx, ny, dLx, dLy)

% ポアソン方程式の係数定義
ac = zeros(nx, ny);
an = zeros(nx, ny);
ae = zeros(nx, ny);
as = zeros(nx, ny);
aw = zeros(nx, ny);

ac(:, :) = 2 * ((1/(dLx)^2) +(1/(dLy)^2));
an(:, :) = 1 / (dLy)^2;
ae(:, :) = 1 / (dLx)^2;
as(:, :) = 1 / (dLy)^2;
aw(:, :) = 1 / (dLx)^2;

% SOR法
eps = 10^(- 14);
ite_max = nx * ny * 20;% 反復回数
alpha = 1.7;% 緩和係数

for ite = 1 : ite_max% SOR法により圧力補正値を求める。

    %ite
    error = 0;

    for i = 1 : nx
        for j = 1 : ny

            if  i == 1% 西側境界条件
                T2d(i, j) = 1;
            elseif  i == nx% 東側境界条件
                T2d(i, j) = T2d(nx - 1, j);
            elseif  j == 1% 南側境界条件
                T2d(i, j) = 0;
            elseif  j == ny% 北側境界条件
                T2d(i, j) = 0;
            else% 内部領域
                T2d_new= ( 1 / ac(i, j) ) * ...
                    ( ae(i, j) * T2d(i + 1, j) + aw(i, j) * T2d(i - 1, j) + an(i, j)* T2d(i, j + 1) + as(i, j)* T2d(i, j - 1) + q2d(i, j)/k );
                error = max(abs(T2d_new - T2d(i, j)), error);
                T2d(i, j) = T2d(i, j) + alpha * (T2d_new - T2d(i, j));
            end
        end
    end

    if error < eps % 収束条件が満たされたらループを抜ける。
        break
    end

    if ite >= ite_max
        disp('最大反復回数に達しました。収束条件を満たしていません。');
    end

end

end

function[u] =  TDMA(a, b, c, d)

%a,cはn-1個の要素
%b,dはn個の要素

n = size(d, 2);% 未知変数の数
e = zeros(n - 1);
f = zeros(n - 1);

% 1番目の係数を求める
e(1) = c(1) / b(1);
f(1) = d(1) / b(1);

% 2からn-1番目の係数を求める
for i = 2 : n - 1
    e(i) = c(i) / (b(i) - a(i - 1) * e(i - 1));
    f(i) = (d(i) - a(i - 1) * f(i - 1)) / (b(i) - a(i - 1) * e(i - 1));
end

% n番目の解を求める
u(n) = (d(n) - a(n - 1) * f(n - 1)) / (b(n) - a(n - 1) * e(n - 1));

% n-1から1番目の解を求める
for i = n - 1 : - 1 : 1
    u(i) = f(i) - e(i) * u(i + 1);
end

end

function[] = vis_contour(T2d, x, y ,Lx, Ly, fignum)

figure(fignum);
h = imagesc(x, y, flip(T2d));
xlabel('y')
ylabel('x')
view(270, 90); % 視点の設定
xticks(0 : round(Lx / 5, 1) : Lx);
yticks(0 : round(Ly / 5, 1) : Ly);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
%axis equal;
axis tight;
axis on;
pbaspect([Lx Ly 1])
colorbar;

end