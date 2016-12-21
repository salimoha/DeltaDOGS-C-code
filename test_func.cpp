/* Calculate the test function at the given point*/
void test_func(point &x)
{
    /* the parabola test function */
    float x0[d], x00=0.3;
    for (int i=0; i<d; i++) x0[i]=x00;
    x.y=0;
    for (int i=0; i<d; i++) x.y+=pow((x.x[i]-x0[i]),2);
}

