function my_array = str2array(my_string)
my_array = strsplit(my_string, {', ',',',' '});
my_array = str2double(my_array);