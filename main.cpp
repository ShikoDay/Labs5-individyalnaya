#include <iostream>
#include <vector>
#include <fstream>

// Класс для решения одномерного уравнения теплопроводности методом прогонки
class HeatConduction1D {
private:
    int N;                   // Количество узлов
    double L;                // Толщина пластины
    double lambda;           // Коэффициент теплопроводности
    double rho;              // Плотность материала
    double c;                // Удельная теплоемкость
    double T0;               // Начальная температура
    double Tl;               // Температура на левом краю
    double Tr;               // Температура на правом краю
    double t_end;            // Время, до которого нужно считать
    double h;                // Шаг по пространственной координате
    double tau;              // Шаг по времени
    std::vector<double> T;   // Вектор для хранения температур
    std::vector<double> alpha, beta;  // Прогоночные коэффициенты

public:
    // Конструктор инициализации параметров задачи
    HeatConduction1D(int nodes, double length, double conductivity, double density, double heatCapacity,
                     double initialTemp, double leftTemp, double rightTemp, double endTime)
        : N(nodes), L(length), lambda(conductivity), rho(density), c(heatCapacity),
          T0(initialTemp), Tl(leftTemp), Tr(rightTemp), t_end(endTime) {
        
        h = L / (N - 1);            // Расчет шага по пространству
        tau = t_end / 100.0;        // Задаем шаг по времени
        T.resize(N, T0);            // Инициализируем вектор температур начальными значениями
        alpha.resize(N, 0.0);       // Вектор коэффициентов alpha для метода прогонки
        beta.resize(N, 0.0);        // Вектор коэффициентов beta для метода прогонки
    }

    // Метод для выполнения численного решения уравнения теплопроводности
    void solve() {  
        double time = 0.0;
        while (time < t_end) {
            time += tau;
            
            // Устанавливаем граничные условия
            alpha[0] = 0.0;
            beta[0] = Tl;

            // Вычисляем прогоночные коэффициенты alpha и beta для внутренней области
            for (int i = 1; i < N - 1; ++i) {
                double ai = lambda / (h * h);
                double bi = 2.0 * lambda / (h * h) + rho * c / tau;
                double ci = lambda / (h * h);
                double fi = -rho * c * T[i] / tau;

                alpha[i] = ai / (bi - ci * alpha[i - 1]);
                beta[i] = (ci * beta[i - 1] - fi) / (bi - ci * alpha[i - 1]);
            }

            // Применяем правое граничное условие
            T[N - 1] = Tr;

            // Обратная прогонка для определения температуры
            for (int i = N - 2; i >= 0; --i) {
                T[i] = alpha[i] * T[i + 1] + beta[i];
            }
        }
    }

    // Метод для записи результатов в файл
    void saveResults(const std::string &filename) const {
        std::ofstream file(filename);
        if (file.is_open()) {
            file << "Распределение температуры на момент времени t = " << t_end << " секунд:\n";
            for (int i = 0; i < N; ++i) {
                file << "x = " << i * h << ", T = " << T[i] << " C\n";
            }
            file.close();
            std::cout << "Результаты сохранены в файл " << filename << std::endl;
        } else {
            std::cerr << "Ошибка открытия файла для записи!" << std::endl;
        }
    }
};

// Основная функция
int main() {
    // Ввод параметров задачи
    int N;
    double L, lambda, rho, c, T0, Tl, Tr, t_end;

    std::cout << "Введите количество узлов по пространственной координате, N: ";
    std::cin >> N;
    std::cout << "Введите конечное время, t_end (с): ";
    std::cin >> t_end;
    std::cout << "Введите толщину пластины, L (м): ";
    std::cin >> L;
    std::cout << "Введите коэффициент теплопроводности, lambda (Вт/(м*С)): ";
    std::cin >> lambda;
    std::cout << "Введите плотность, rho (кг/м^3): ";
    std::cin >> rho;
    std::cout << "Введите удельную теплоемкость, c (Дж/(кг*С)): ";
    std::cin >> c;
    std::cout << "Введите начальную температуру, T0 (С): ";
    std::cin >> T0;
    std::cout << "Введите температуру на левом краю, Tl (С): ";
    std::cin >> Tl;
    std::cout << "Введите температуру на правом краю, Tr (С): ";
    std::cin >> Tr;

    // Создаем объект и запускаем решение
    HeatConduction1D heatConduction(N, L, lambda, rho, c, T0, Tl, Tr, t_end);
    heatConduction.solve();

    // Сохраняем результаты
    heatConduction.saveResults("temperature_results.txt");

    return 0;
}
