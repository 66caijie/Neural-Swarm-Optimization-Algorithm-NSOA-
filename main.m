%_________________________________________________________________________%
%                                                                         %
%  Neural Swarm Optimization Algorithm (NSOA) source codes version 1.0    %
%                                                                         %
%  Developed in:	MATLAB 9.5 (R2018b)                                   %
%                                                                         %
%  Programmer:		Jie Cai                                               %
%                                                                         %
%  Original paper:	Jie Cai, Jiayao Chen                                  %
%					A novel metaheuristic algorithm for solving           %
%                   global optimization and engineering problems.         %
%                                                                         %
%_________________________________________________________________________%
clc
clear all
close all

N = 60; % Number of search agents
Max_iter = 1000; % Maximum number of iterations
F_name = 'F'; % Name of the test function (ignored in Get_Functions_details)

% Load details of the selected benchmark function
[lb, ub, dim, fobj] = Get_Functions_details(F_name);

tic;
[Fbest, Xbest, CNVG] = NSOA(N, Max_iter, lb, ub, dim, fobj);
Run_time = toc;

semilogy(1:Max_iter, CNVG, 'color', 'r', 'linewidth', 2.5);
title('Convergence curve');
xlabel('Iteration');
ylabel('Best score obtained so far');
display(['The running time is: ', num2str(Run_time)]);
display(['The best fitness is: ', num2str(Fbest)]);
display(['The best position is: ', num2str(Xbest)]);
