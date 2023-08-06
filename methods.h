#pragma once

//прототипы
void Adams_J(vector<double>* J, double x, double y, unsigned long long* count, double t);
void Euler_J(vector<double>* J, double x, double y, unsigned long long* count, double t);
double R(double x, double y, unsigned long long* counter);
double main_diag(double x, double y, unsigned long long* counter);
double side_diag(double x, double y, unsigned long long* counter);
void Newton_method(int i, vector<point>* P, double (*F[])(int, vector<point>*, ARG, double), double* X,
	point left_border, point right_border, double acc, unsigned long long* count, bool flag, double tau);

//функции системы
double R(double x, double y, unsigned long long* counter) {
	*counter += 6;
	return x * (x - 2.0) * (6.0 + 2.0 * x * (y - 1.0) - y * (10.0 - 3 * y));
};
//элемент главной диагонали для матрицы Якоби
double main_diag(double x, double y, unsigned long long* counter) {
	*counter += 7;
	return 6.0 - 3.0 * x * x * (y - 1.0) + y * (3.0 * y - 10.0) + x * (y * (14.0 - 3.0 * y) - 10.0);
};
//элемент побочной диагонали
double side_diag(double x, double y, unsigned long long* counter) {
	*counter += 3;
	return x * (2.0 - x) * (x + 3.0 * y - 5.0);
};

//метод Ньютона для системы из двух уравнений (допустимая область задается прямоугольником, начальная точка должна в него попадать)
//аргументы:
//индекс, массив точек, массив функций, обратная матрица Якоби (развернутая в строку), 
//левая нижняя вершина допустимой области, правая верхняя вершина допустимой области, точность, счетчик умножений, флаг (система Эйлера - 0 или Адамса - 1)
void Newton_method(int i, vector<point>* P, double (*F[])(int, vector<point>*, ARG, double), double* X, 
	point left_border, point right_border, double acc, unsigned long long* count, bool flag, double tau) {

	//промежуточные точки
	double Temp[2]{ (*P)[i - 1].x, (*P)[i - 1].y };
	double NewTemp[2];
	double TP[2];
	vector<point> trac;
	//инициализация начальной точки (в данном массиве временная коорд 1 - точка внутри допустимой области,0 - вне допусти. области
	trac.push_back({ (*P)[i - 1].x, (*P)[i - 1].y, 1 });
	//индекс последней точки в траектории, принадлежащей допустимой области
	unsigned int trac_index = 0;
	//флаги выхода за границу по оси x и y
	bool cross_x, cross_y;
	//параметр
	double px = 1.0;
	double py = 1.0;

	//расстояние между двумя точками в евклидовой норме
	auto norma = [](double X[2], double Y[2], unsigned long long* count) {
		double x = X[0] - Y[0];
		double y = X[1] - Y[1];
		*count += 2;
		return sqrt(x * x + y * y);
	};

	int counter = 0;

	//обратная матрица Якоби
	vector<double> J;
	J.push_back(1);
	J.push_back(0);
	J.push_back(0);
	J.push_back(1);

	//цикл с пост условием (пока не выполнено условие по норме)
	do {
		if (flag)
			Adams_J(&J, Temp[0], Temp[1], count, tau);
		else
			Euler_J(&J, Temp[0], Temp[1], count, tau);

		//расчет новой точки
		NewTemp[0] = Temp[0] - px * ((*F[0])(i, P, Temp[0], Temp[1], count, tau) * J[0] +
			(*F[1])(i, P, Temp[0], Temp[1], count, tau) * J[1]);
		NewTemp[1] = Temp[1] - py * ((*F[0])(i, P, Temp[0], Temp[1], count, tau) * J[2] +
			(*F[1])(i, P, Temp[0], Temp[1], count, tau) * J[3]);
		*count += 6;



		//отмечаем точку в траектории
		trac.push_back({ NewTemp[0], NewTemp[1],0 });

		//проверка выхода за границы
		cross_x = (NewTemp[0] > right_border.x) || (NewTemp[0] < left_border.x);
		cross_y = (NewTemp[1] > right_border.y) || (NewTemp[1] < left_border.y);

		if (cross_x && cross_y) {
			px = -0.9 * px;
			py = -0.9 * py;
			*count += 2;
			Temp[0] = trac[trac_index].x;
			Temp[1] = trac[trac_index].y;
			continue;
		}
		else if (!cross_x && cross_y) {
			py = -0.9 * py;
			*count += 1;
			Temp[0] = trac[trac_index].x;
			Temp[1] = trac[trac_index].y;
			continue;
		}
		else if (cross_x && !cross_y) {
			px = -0.9 * px;
			*count += 1;
			Temp[0] = trac[trac_index].x;
			Temp[1] = trac[trac_index].y;
			continue;
		}
		else if (!cross_x && !cross_y) {
			trac_index = trac.size() - 1;
			trac[trac_index].t = 1;
		};

		//слишком малые py и px нежелательны
		if (abs(py) < 0.01) py = py * 100;
		if (abs(px) < 0.01) px = px * 100;

		//cout << p << endl;
		//обновление точки
		TP[0] = Temp[0];
		TP[1] = Temp[1];
		Temp[0] = NewTemp[0];
		Temp[1] = NewTemp[1];
		counter++;

		//условие прекращения цикла
		//значения компонент вектора функций меньше eps
		//расстояние между точками меньше eps
	} while (
		((abs((*F[0])(i, P, NewTemp[0], NewTemp[1], count, tau)) > acc) || (abs((*F[1])(i, P, NewTemp[0], NewTemp[1], count, tau)) > acc)) &&
		(norma(TP, NewTemp, count) > acc)
		);

	//запись результата
	X[0] = NewTemp[0];
	X[1] = NewTemp[1];

	//for (int i = 0; i < trac.size(); i++) {	cout << "point" << i << " = (" << trac[i].x << ", " << trac[i].y << ");  " << trac[i].t << endl;};
};

//обратная матрица Якоби для системы Эйлера
void Euler_J(vector<double>* J, double x, double y, unsigned long long* count, double t) {

	*count += 11;

	double md = main_diag(x, y, count);
	double frac = 2 * t;

	double a = 1 + frac * md;
	double b = frac * side_diag(x, y, count);
	double c = -frac * side_diag(y, x, count);
	double d = 1 - frac * md;

	double det = a * d - b * c;

	//результат
	(*J)[0] = d / det;
	(*J)[1] = -b / det;
	(*J)[2] = -c / det;
	(*J)[3] = a / det;

};

//обратная матрица Якоби для системы Адамса (4 пор)
void Adams_J(vector<double>* J, double x, double y, unsigned long long* count, double t) {

	*count += 12;

	double frac = t * 251.0 / 360.0;
	double md = main_diag(x, y, count);
	double a = 1 + frac * md;
	double b = frac * side_diag(x, y, count);
	double c = -frac * side_diag(y, x, count);
	double d = 1 + frac * md;

	double det = a * d - b * c;

	//результат
	(*J)[0] = d / det;
	(*J)[1] = -b / det;
	(*J)[2] = -c / det;
	(*J)[3] = a / det;
};