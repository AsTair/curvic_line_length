#include <cmath>
#include <algorithm>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <string>

using namespace std;

#define ARG double,double,unsigned long long*

//тип точки
struct point {
	//координата абцисс
	double x;
	//координата ординат
	double y;
	//временная координата
	double t;
};

//пользовательские функции
#include"methods.h"
#include"get_trajectory.h"

int main(){
	
	//глобальный счетчик операций умножения
	unsigned long long mult_counter[6]{ 0,0,0,0,0,0 };

	//массивы точек
	vector<point> points1;
	points1.push_back({ 1, 1, 0.0 });

	vector<point> points2;
	points2.push_back({ 1, 1, 0.0 });
	
	vector<point> points3;
	points3.push_back({ 1, 1, 0.0 });
	
	vector<point> points4;
	points4.push_back({ 1, 1, 0.0 });
	
	vector<point> points5;
	points5.push_back({ 1, 1, 0.0 });

	vector<point> points6;
	points6.push_back({ 1, 1, 0.0 });

	//расстояние между двумя точками в евклидовой норме
	auto norma = [](point A, point B, unsigned long long* count) {
		double x = A.x - B.x;
		double y = A.y - B.y;
		*count += 2;
		return sqrt(x * x + y * y);
	};
	
	//для таймера
	time_t start, end;
	double sec[6]{0,0,0,0,0,0};
	double t[6]{ 1E-5, 2E-5, 1E-6, 2E-6, 1E-7, 2E-7 }; double step_number = 6;
	double eps[6]{ 1E-9, 1E-10, 1E-11, 1E-12, 1E-13, 1E-14 }; double eps_number = 6;
	vector<point>* trac[6]{ &points1, &points2, &points3, &points4, &points5, &points6 };

	//макс. кол-во точек
	const unsigned long long N = 200000000;
	unsigned long long temp;
	
	//вывод данных
	ofstream data;
	
	//расчет длины
	double length[6]{ 0,0,0,0,0,0 };

	//погрешность решения ОДУ (правило Рунге)
	double acc[3]{ 0,0,0 };
	
	for (int j = 0; j < eps_number; j++) { //для разных точностей
		
		cout << "Accuracy = " << eps[j] << endl;
		
		for (int i = 0; i < 6; i++) { //для шести разных шагов
			cout << "start " << i + 1 << " trajectory" << endl;
			time(&start);
			get_trajectory(trac[i], t[i], 2*t[i], eps[j], N, &mult_counter[i]);
			time(&end);
			sec[i] = difftime(end, start);
			cout << "finish 1 trajectory, time = " << sec[i] << endl << endl;
		};
				
		//вычисление длины
		for (int k = 0; k < step_number; k++) {
			length[k] = 0;

			for (int i = 1; i < (*trac[k]).size(); i++) {

				length[k] += norma((*trac[k])[i], (*trac[k])[i - 1], &mult_counter[k]);

			};
			length[k] += norma((*trac[k])[0], (*trac[k])[(*trac[k]).size() - 1], &mult_counter[k]);
		};

		//вычисление погрешности ОДУ
		acc[0] = norma((*trac[0])[(*trac[0]).size() - 1], (*trac[1])[(*trac[1]).size() - 1], &temp);
		acc[1] = norma((*trac[2])[(*trac[2]).size() - 1], (*trac[3])[(*trac[3]).size() - 1], &temp);
		acc[2] = norma((*trac[4])[(*trac[4]).size() - 1], (*trac[5])[(*trac[5]).size() - 1], &temp);

		//запись основных данных
		data.open("data_" + to_string(j+1) + ".txt");

		data << "Accuracy1 = " << acc[0]/15.0 << endl;
		data << "Accuracy2 = " << acc[1]/15.0 << endl;
		data << "Accuracy3 = " << acc[2]/15.0 << endl << endl << endl;
		
		data << "eps = " << eps[j] << endl << endl;;

		for (int i = 0; i < 6; i++) {
			data << "#" << i + 1 << endl;
			data << t[i] << "   &   " << (*trac[i]).size() << "   &   " << mult_counter[i] << "   &   " << length[i] << "  \\ " << endl;
			data << "time = " << sec[i] << " s \n";
		};

		data.close();
		/*
		//запись траектории
		for (int i = 0; i < 2; i++) {
			data.open("DATA_" + to_string(j) + "_" + to_string(i + 1) + ".dat");
			data << "{";
			for (int j = 0; j < (*trac[i]).size() - 1; j++) {
				data << "{" << (*trac[i])[j].x << "," << (*trac[i])[j].y << "},";
			};
			data << "{" << (*trac[i])[(*trac[i]).size() - 1].x << "," << (*trac[i])[(*trac[i]).size() - 1].y << "}} \n";
			data.close();
		};*/

		//очистка траекторий
		for (int i = 0; i < 6; i++) {
			(*trac[i]).clear();
			(*trac[i]).shrink_to_fit();
			(*trac[i]).push_back({ 1,1,0.0 });
		};

	};

	return 0;
}

