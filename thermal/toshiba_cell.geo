Merge "toshiba_cell.stp";
//+
Physical Surface("terminal_negative", 49) = {11};
//+
Physical Surface("terminal_positive", 50) = {12};
//+
Physical Surface("housing_I", 51) = {1};
//+
Physical Surface("housing_VI", 52) = {6};
//+
Physical Surface("housing_II", 53) = {2};
//+
Physical Surface("housing_V", 54) = {4};
//+
Physical Surface("housing_III", 55) = {5};
//+
Physical Surface("housing_IV", 56) = {3};
//+
Physical Surface("negleted", 57) = {7, 10, 9, 8, 15, 16, 13, 14};

Surface Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
Volume(1) = {1};//+
Physical Volume("volume", 58) = {1};
