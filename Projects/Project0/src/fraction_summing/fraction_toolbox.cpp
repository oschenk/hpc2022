#include <iostream>

#include "fraction_toolbox.hpp"

void print_fraction(fraction frac)
{
    std::cout << frac.num << '/' << frac.denom << std::endl;
}

void print_fraction_array(fraction frac_array[], int n)
{
    std::cout << "[ " << frac_array[0].num << '/' << frac_array[0].denom << std::endl;
    for (int i = 1; i < n-1; i++)
    {
        std::cout << "  ";
        print_fraction(frac_array[i]);
    }
    std::cout << "  " << frac_array[n-1].num << '/' << frac_array[n-1].denom << " ]" << std::endl;
}

fraction square_fraction(fraction frac)
{
    //TODO: implement function 2
}

//TODO: implement function 3


double fraction2double(fraction frac)
{
    //TODO: implement function 4
}

int gcd(int a, int b)
{
    //TODO: implement function 5
}

//TODO: implement function 6


void reduce_fraction_inplace(fraction & frac)
{
    //TODO: implement function 7

    //TODO: add short comment to explain which of the gcd() functions your code is calling
}

fraction add_fractions(fraction frac1, fraction frac2)
{
    //TODO: implement function 8
}

double sum_fraction_array_approx(fraction frac_array[], int n)
{
    //TODO: implement function 9
    
    //TODO: add short comment to explain why this function is approximate
}

//TODO: implement function 10


void fill_fraction_array(fraction frac_array[], int n)
{
    fraction temp_frac;
    temp_frac.num = 1;
    for (int i = 1; i <= n; i++)
    {
        temp_frac.denom = i * (i+1);
        frac_array[i-1] = temp_frac;
    }
}

