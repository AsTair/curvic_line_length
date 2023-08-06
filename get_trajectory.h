#pragma once

void get_trajectory(vector<point>* points, double tau, double h, double eps, int N, unsigned long long* mult_counter ){

	double Xk[2];

	//рассто€ние между двум€ точками в евклидовой норме
	auto norma = [](point A, point B, unsigned long long* count) {
		double x = A.x - B.x;
		double y = A.y - B.y;
		*count += 2;
		return sqrt(x * x + y * y);
	};

/*	//не€вный метод Ёйлера
	//л€мбда фукнции дл€ не€вного метода Ёйлера
	//задание функций системы
	auto euler_f1 = [](int i, vector<point>* P, double x, double y, unsigned long long* count, double t) {
		*count += 1;
		return x - (*P)[i - 1].x - t * R(x, y, count);
	};
	auto euler_f2 = [](int i, vector<point>* P, double x, double y, unsigned long long* count, double t) {
		*count += 1;
		return y - (*P)[i - 1].y + t * R(y, x, count);
	};

	//массив фукнций
	double (*euler_F[])(int, vector<point>*, ARG, double) = { euler_f1, euler_f2 };

	//расчет X1,X2,X3 по не€вному методу Ёйлера
	for (int i = 1; i < 4; i++) {
		Newton_method(i, points, euler_F, Xk, { (*points)[i - 1].x - h,(*points)[i - 1].y - h,0 }, { (*points)[i - 1].x + h,(*points)[i - 1].y + h,0 }, eps, mult_counter, 0, tau);
		(*points).push_back({ Xk[0],Xk[1],i * tau });
	}; */

	//метод –унге- утты 4 пор€дка дл€ системы из двух уравнений
	double k1, k2, k3, k4; //дл€ х
	double q1, q2, q3, q4; //дл€ у
	for (int i = 1; i < 4; i++) {

		k1 = tau * R((*points)[i - 1].x, (*points)[i - 1].y, mult_counter);
		q1 = -tau * R((*points)[i - 1].y, (*points)[i - 1].x, mult_counter);

		k2 = tau * R((*points)[i - 1].x + 0.5 * k1, (*points)[i - 1].y + 0.5 * q1, mult_counter);
		q2 = -tau * R((*points)[i - 1].y + 0.5 * q1, (*points)[i - 1].x + 0.5 * k1, mult_counter);

		k3 = tau * R((*points)[i - 1].x + 0.5 * k2, (*points)[i - 1].y + 0.5 * q2, mult_counter);
		q3 = -tau * R((*points)[i - 1].y + 0.5 * q2, (*points)[i - 1].x + 0.5 * k2, mult_counter);

		k4 = tau * R((*points)[i - 1].x + k3, (*points)[i - 1].y + q3, mult_counter);
		q4 = -tau * R((*points)[i - 1].y + q3, (*points)[i - 1].x + k3, mult_counter);

		(*points).push_back({ (*points)[i - 1].x + (k1 + k4 + 2.0 * (k2 + k3)) / 6.0,
			(*points)[i - 1].y + (q1 + q4 + 2.0 * (q2 + q3)) / 6.0,
			tau * i });

		*mult_counter += 21;
	};

	//не€вный метод јдамса 4 пор€дка
	auto adams_f1 = [](int i, vector<point>* P, double x, double y, unsigned long long* count, double t) {

		*count += 7;

		double R0 = R((*P)[i - 4].x, (*P)[i - 4].y, count);
		double R1 = R((*P)[i - 3].x, (*P)[i - 3].y, count);
		double R2 = R((*P)[i - 2].x, (*P)[i - 2].y, count);
		double R3 = R((*P)[i - 1].x, (*P)[i - 1].y, count);
		double R4 = R(x, y, count);

		double RR = t * (251.0 * R4 + 646.0 * R3 - 264 * R2 + 106 * R1 - 19 * R0) / 720.0;

		return x - (*P)[i - 1].x - RR;
	};
	auto adams_f2 = [](int i, vector<point>* P, double x, double y, unsigned long long* count, double t) {

		*count += 7;

		double R0 = R((*P)[i - 4].y, (*P)[i - 4].x, count);
		double R1 = R((*P)[i - 3].y, (*P)[i - 3].x, count);
		double R2 = R((*P)[i - 2].y, (*P)[i - 2].x, count);
		double R3 = R((*P)[i - 1].y, (*P)[i - 1].x, count);
		double R4 = R(y, x, count);

		double RR = t * (251.0 * R4 + 646.0 * R3 - 264 * R2 + 106 * R1 - 19 * R0) / 720.0;

		return y - (*P)[i - 1].y + RR;
	};
	//массив фукнций
	double (*adams_F[])(int, vector<point>*, ARG, double) = { adams_f1, adams_f2 };

	//расчет по не€вному методу јдамса
	int i = 4;
	do {
		Newton_method(i, points, adams_F, Xk, { (*points)[i - 1].x - h,(*points)[i - 1].y - h,0 }, { (*points)[i - 1].x + h,(*points)[i - 1].y + h,0 }, eps, mult_counter, 1, tau);
		(*points).push_back({ Xk[0],Xk[1],i * tau });
		*mult_counter += 1;
		i++;
	} while (((*points).size() < N) && (norma((*points)[0], (*points)[i - 1], mult_counter) > h));

};

