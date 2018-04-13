namespace PSE{
    unsigned factorial(unsigned n)
    {
        if (n == 1)
            return 1;
        else
            return n * factorial(n - 1);
    }

}
