#!/bin/bash
Nt=( 10 20 50 100)
omega=( 5.0 2.5 1.0 0.5)
for ((i=0;i<4;i+=1))
do
	./a.out ${Nt[i]} ${omega[i]}
done
