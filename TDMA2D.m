% ������
clear all;

% �p�����[�^�[
Lx = 1;% �v�Z�̈�
Ly = 1;% �v�Z�̈�
nx = 20;% �v�f��
ny = 20;% �v�f��
dLx = Lx / nx;% �v�f�T�C�Y
dLy = Ly / ny;% �v�f�T�C�Y
k = 1;% �M�`����

% ���W�̐���
x = 0 : Lx / nx : Lx;
y = Ly : - Ly / ny : 0;
[X, Y] = meshgrid(x, y);

% �z��m��
T2d = zeros(nx, ny);

% �M���̐ݒ�
q2d = zeros(nx, ny);
%q2d(round(nx / 2), round(ny / 2)) = 1000;% �����ɒP�ʔM����ݒu����

% ���Ԍv���J�n
tic

% �A���������̉�@��I��
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

        % �����@�̃p�����[�^�[
        eps = 10^(-14);% ���e�덷
        ite_max = nx * ny *20;% �ő唽����
        omega = 1.3;% �ɘa�W��

        for ite = 1 : ite_max

            error = 0;

            %ite

            for i = 2 : nx - 1

                % �}x�����E����
                T2d(1, :) = 1;% �f�B���N������
                T2d(end, :) = T2d(end - 1, :) ; % �m�C�}������
                %T2d(end, :) = 0; % �f�B���N������

                % j = 1 �̎� (-y���̋��E����)
                b(1) = 1;
                c(1) = 0;% �f�B���N������
                % c(1) = -1;% �m�C�}������
                d(1) = 0;

                % j = 2 ~ ny - 1 �̎�(�����̈�)
                for j = 2 : ny - 1
                    a(j - 1) = beta;
                    b(j) = gamma / omega;
                    c(j) = beta;
                    d(j) = - alpha * T2d(i + 1, j) - alpha * T2d(i - 1, j) - q2d(i, j) / k + (1 - omega) * gamma * T2d(i, j)/ omega;
                end

                % j = ny �̎�(+y���̋��E����)
                a(ny - 1) = 0;% �f�B���N������
                % a(ny - 1) = -1;% �m�C�}������
                b(ny) = 1;
                d(ny) = 0;

                % TDMA
                T2d_new(i, :) = TDMA(a, b, c, d);

                % �덷�̌v�Z
                for j = 2 : ny - 1
                    cor = T2d(i, j) - T2d_new(i, j);
                    error = max(error, abs(cor));
                    T2d(i, j) = T2d_new(i, j);
                end

                % �p�����͋��E�����̕��ςƂ���
                T2d(1, 1) = (T2d(1, 2) + T2d(2, 1)) / 2;
                T2d(1, end) = (T2d(1, end - 1) + T2d(2, end)) / 2;
                T2d(end, 1) = (T2d(end, 2) + T2d(end - 1, 1)) / 2;
                T2d(end, end) = (T2d(end - 1, end) + T2d(end, end - 1)) / 2;

            end

            if error < eps %&& ite> 100 % �����������������ꂽ�烋�[�v�𔲂���B
                break
            end

            if ite >= ite_max
                disp('�ő唽���񐔂ɒB���܂����B���������𖞂����Ă��܂���B');
            end

        end

    otherwise

        error("�I�����ꂽ��@�F" + method + "�̓T�|�[�g���Ă��Ȃ��ł��B" + ...
            "SOR,TDMA����I�����Ă��������B");

end

% ���Ԍv���I��
toc

% ����
vis_contour(T2d, x, y, Lx, Ly, 2)

%% �ȉ��֐�
function[T2d] = SOR(T2d, q2d, k, nx, ny, dLx, dLy)

% �|�A�\���������̌W����`
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

% SOR�@
eps = 10^(- 14);
ite_max = nx * ny * 20;% ������
alpha = 1.7;% �ɘa�W��

for ite = 1 : ite_max% SOR�@�ɂ�舳�͕␳�l�����߂�B

    %ite
    error = 0;

    for i = 1 : nx
        for j = 1 : ny

            if  i == 1% �������E����
                T2d(i, j) = 1;
            elseif  i == nx% �������E����
                T2d(i, j) = T2d(nx - 1, j);
            elseif  j == 1% �쑤���E����
                T2d(i, j) = 0;
            elseif  j == ny% �k�����E����
                T2d(i, j) = 0;
            else% �����̈�
                T2d_new= ( 1 / ac(i, j) ) * ...
                    ( ae(i, j) * T2d(i + 1, j) + aw(i, j) * T2d(i - 1, j) + an(i, j)* T2d(i, j + 1) + as(i, j)* T2d(i, j - 1) + q2d(i, j)/k );
                error = max(abs(T2d_new - T2d(i, j)), error);
                T2d(i, j) = T2d(i, j) + alpha * (T2d_new - T2d(i, j));
            end
        end
    end

    if error < eps % �����������������ꂽ�烋�[�v�𔲂���B
        break
    end

    if ite >= ite_max
        disp('�ő唽���񐔂ɒB���܂����B���������𖞂����Ă��܂���B');
    end

end

end

function[u] =  TDMA(a, b, c, d)

%a,c��n-1�̗v�f
%b,d��n�̗v�f

n = size(d, 2);% ���m�ϐ��̐�
e = zeros(n - 1);
f = zeros(n - 1);

% 1�Ԗڂ̌W�������߂�
e(1) = c(1) / b(1);
f(1) = d(1) / b(1);

% 2����n-1�Ԗڂ̌W�������߂�
for i = 2 : n - 1
    e(i) = c(i) / (b(i) - a(i - 1) * e(i - 1));
    f(i) = (d(i) - a(i - 1) * f(i - 1)) / (b(i) - a(i - 1) * e(i - 1));
end

% n�Ԗڂ̉������߂�
u(n) = (d(n) - a(n - 1) * f(n - 1)) / (b(n) - a(n - 1) * e(n - 1));

% n-1����1�Ԗڂ̉������߂�
for i = n - 1 : - 1 : 1
    u(i) = f(i) - e(i) * u(i + 1);
end

end

function[] = vis_contour(T2d, x, y ,Lx, Ly, fignum)

figure(fignum);
h = imagesc(x, y, flip(T2d));
xlabel('y')
ylabel('x')
view(270, 90); % ���_�̐ݒ�
xticks(0 : round(Lx / 5, 1) : Lx);
yticks(0 : round(Ly / 5, 1) : Ly);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
%axis equal;
axis tight;
axis on;
pbaspect([Lx Ly 1])
colorbar;

end