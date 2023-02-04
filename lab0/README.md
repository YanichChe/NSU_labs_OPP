__Практическое задание №1 “Скалярное произведение векторов”__

I Написать 3 программы, каждая из которых вычисляет скалярное произведение двух векторов:

a) последовательная программа

b) параллельная, использующая коммуникации типа точка-точка (MPI_Send, MPI_Recv)

c) параллельная, использующая коллективные коммуникации (MPI_Scatter, MPI_Reduce)

II Замерить время работы последовательной программы и параллельных на 2, 4, 8, 16, 24 процессах. Рекомендуется провести несколько замеров для каждого варианта запуска и выбрать минимальное время.

III Построить графики времени, ускорения и эффективности.

IV Составить отчет, содержащий исходные коды разработанных программ и построенные графики.

_Требования:_
- длину векторов выбирать таким образом, чтобы время работы последовательной программы было не менее 30 сек;
- в параллельных программах изначально векторы должны полностью инициализироваться на 0-м процессе. Для параллельного расчета 0-й процесс должен раздавать части векторов остальным;
- в параллельных программах полное скалярное произведение должно в результате выводиться на экран 0-м процессом.