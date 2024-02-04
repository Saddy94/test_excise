## A second-level heading
Задача

Нужно написать программу на python, которая рассчитывает кватернион целевой ориентации КА на низкой околоземной орбите. Целевая ориентация состоит в наведении одной из осей связанной системы координат (ССК) КА на заданную точку на поверхности Земли (подробнее ниже).
Требования к программе
В программе должна быть возможность указать интересующий промежуток времени.
Программа должна формировать кватернион целевой ориентации для каждого момента времени внутри заданного промежутка с задаваемым шагом (по умолчанию 1 секунда).
Эти и другие необходимые параметры должны задаваться хардкодом в начале главного скрипта программы.
Программа не должна запрашивать никакого ввода от пользователя и должна работать без его вмешательства от запуска главного скрипта до завершения.
Результаты работы программа должна записать в файл в описанном ниже формате.
## A second-level heading
Системы координат
## A third-level heading
ССК
Целевая ориентация ССК описывается следующий образом:
Ось Y направлена от целевой точки к аппарату.
Ось Z лежит в плоскости, образуемой нормалью к орбите и осью Y, и направлена в сторону от нормали к орбите.
Ось X дополняет тройку до правой и направлена в сторону скорости аппарата.
Кватернион целевой ориентации КА (обозначается как q) должен описывать ориентацию базиса ССК относительно базиса инерциальной системы координат (ИСК) так, что проекции произвольного вектора v на оси ИСК (обозначается как vi) выражаются через проекции этого же вектора на оси ССК (обозначается как vs) следующим образом:
vi = q ∘ vs ∘  ͞q,
где чертой обозначается сопряжённый кватернион.
### A third-level heading
ИСК
Инерциальной системой координат является классическая ECI J2000, которую можно считать равной системе ICRF.
Формат файла с результатами
Файл должен содержать по одной строке на один момент времени. В каждой строке должно быть 11 чисел (дробную часть отделяет точка, не запятая), разделённых одним символом пробела:
Время в формате количества юлинский дней в шкале TT, прошедших от эпохи J2000.
X-компонента радиус-вектора КА в осях ИСК.
Y-компонента радиус-вектора КА в осях ИСК.
Z-компонента радиус-вектора КА в осях ИСК.
X-компонента скорости КА в осях ИСК.
Y-компонента скорости КА в осях ИСК.
Z-компонента скорости КА в осях ИСК.
Скалярная часть кватерниона.
X-компонента векторной части кватерниона.
Y-компонента векторной части кватерниона.
Z-компонента векторной части кватерниона.
Исходные данные
Координаты КА в ИСК на начальный момент времени (м):
(4362521.19692133, -2174459.71448059, 4720847.40402189)

Скорость КА в ИСК на начальный момент времени (м/с):
(5356.39915538069, 4741.41348686709, -2761.5472632395)

Начальный момент времени (шкала TT):
Юлианские дни от эпохи J2000: 8084.05644194318
Или в календарном формате:  2022-02-18 13:21:16.584

Географические координаты целевой точки на Земле:
(45.920266° с.ш., 63.342286° з.д.)
