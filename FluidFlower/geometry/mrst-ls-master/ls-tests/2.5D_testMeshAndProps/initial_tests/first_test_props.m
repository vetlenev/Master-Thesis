%% Test how to import fluid properties from ECLIPSE deck
clc; close all force; clear
mrstModule add ad-props deckformat 

deck = readEclipseDeck('C:\Users\lsalo\matlab\mrst-2019a\myprojects\2.5Dtestmesh\ftest.DATA');
deck = convertDeckUnits(deck);
fluid = initDeckADIFluid(deck);