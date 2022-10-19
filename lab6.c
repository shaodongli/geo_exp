/*
 code in lab 6
 references
 https://tex.stackexchange.com/questions/348651/c-code-to-add-in-the-document
 https://stackoverflow.com/questions/65923375/gnuplot-multiplot-and-insets
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
typedef struct Model Model;

struct Model
{
    double h; // height of a single layer
    double v;
    double d; // depth of the layer
};

int get_horizontal_distance_and_traveltime(double depth, double theta, double v1, double v2, double h1, double * x, double * t)
{
    // this function assumes depth > h1
    *x = (depth - h1) * tan(theta);
    double p = sin(theta) / v2; // ray parameter
    double eta1 = sqrt(1.0 / v1 / v1 - p * p); // vertical slowness in layer 1
    double eta2 = sqrt(1.0 / v2 / v2 - p * p); // vertical slowness in layer 2
    double theta1 = asin(p * v1);
    *x += h1 * tan(theta1);
    *t = (*x) * p + h1 * eta1 + (depth - h1) * eta2;
    return 0;
}

int get_index_in_Model(const double depth, Model * vel, const int n)
{
    int index = -1;
    
    for (int i = 0; i < n - 1; i++)
    {
        if ( (depth - vel[i].d) * (depth - vel[i+1].d) < 0)
        {
            index = i;
            break;
        }
            
    }
    return index;
}

int get_horizontal_distance_and_traveltime_Model(double depth, double theta, Model * vel, const int n, double * x, double * t)
{
    int ni = get_index_in_Model(depth, vel, n);
    if (ni < 0)
        return ni;
    // this function assumes depth > h1
    *x = 0;
    *t = 0;
    double p;
    double eta;
    double eta2;
    double theta2;

    *x += (depth - vel[ni].d) * tan(theta);
    p = sin(theta) / vel[ni].v;
    eta = sqrt(1.0 / vel[ni].v / vel[ni].v - p * p);
    *t += eta * (depth - vel[ni].d);
    
    for (int i = ni; i > 0; i--)
    {
        theta = asin(p * vel[i-1].v);
        eta = sqrt(1.0 / vel[i-1].v / vel[i-1].v - p * p);
        *x += (vel[i].d - vel[i-1].d) * tan(theta);
        *t += eta * (vel[i].d - vel[i-1].d);
    }
    *t += (*x) * p;
    return 0;
}

double get_take_off_angle(double depth, double v1, double v2, double h1, double Delta, double * t)
{
    double x = 1e6;
    double tl = 0.0, tu = 1.5707;
    double epsilon = 1e-5;
    double tm = tu;
    while (fabs(x - Delta) > epsilon)
    {
        tm = (tu + tl) * 0.5;
        get_horizontal_distance_and_traveltime(depth, tm, v1, v2, h1, &x, t);
        if (x > Delta)
            tu = tm;
        else
            tl = tm;
    }
    return tm;
}

double get_take_off_angle_Model(double depth, Model * vel, const int n, double Delta, double * t)
{
    double x = 1e6;
    double tl = 0.0, tu = 1.5707;
    double epsilon = 1e-2;
    double tm = tu;
    while (fabs(x - Delta) > epsilon)
    {
        tm = (tu + tl) * 0.5;
        get_horizontal_distance_and_traveltime_Model(depth, tm, vel, n, &x, t);
        if (x > Delta)
            tu = tm;
        else
            tl = tm;
    }
    return tm;
}

#define PI 3.141592653589793
double rad2deg(const double radian)
{
    return radian / PI * 180.0;
}

int plot_traveltime(Model * vel, const int n, const double depth)
{
    const int num = 1000;
    double x, t;
    double theta;
    double di = PI * 0.499 / 1000.0;
    FILE * fp = fopen("pg.dat", "w");
    if (fp != NULL)
    {
        for (int i = 1; i < num; i++)
        {
            theta = di * i;
            get_horizontal_distance_and_traveltime_Model(depth, theta, vel, n, &x, &t);
            fprintf(fp, "%lf %lf\n", x, t);
            if (x > 999)
                break;
        }
        fclose(fp);
    }
    else
        perror("pg.dat");
    
    FILE * fp0 = fopen("pg.sh", "w");
    if (fp0 != NULL)
    {
        fprintf(fp0, "set term pdf\n");
        fprintf(fp0, "set output \"pg.pdf\"\n");
        fprintf(fp0, "set size ratio 1.\n");
        fprintf(fp0, "%s\n", "set xlabel \"{/Times-Italic {/Symbol D}} ({/Times-Italic km})\"");
        fprintf(fp0, "%s\n", "set ylabel \"{/Times-Italic Traveltime}  ({/Times-Italic s})\"");
        fprintf(fp0, "set xrange [0:1000]\n");
        fprintf(fp0, "plot \"pg.dat\" with lines notitle lc \"black\"\n");
        fclose(fp0);
    }
    else
        perror("pg.sh");
    
    fp0 = fopen("pg2.sh", "w");
    if (fp0 != NULL)
    {
        fprintf(fp0, "set term pdf\n");
        fprintf(fp0, "set output \"pg2.pdf\"\n");
        fprintf(fp0, "set key top left\n");
        fprintf(fp0, "set multiplot layout 1,2\n");
        fprintf(fp0, "s1x = 0.6\n");
        fprintf(fp0, "s1y = 1.0\n");
        fprintf(fp0, "set size s1x, s1y\n");
        fprintf(fp0, "o1x = 0.0\n");
        fprintf(fp0, "o1y = 0.0\n");
        fprintf(fp0, "set origin o1x, o1y\n");
        // subplot 1
        fprintf(fp0, "%s\n", "set xlabel \"{/Times-Italic {/Symbol D}} ({/Times-Italic km})\"");
        fprintf(fp0, "%s\n", "set ylabel \"{/Times-Italic Traveltime}  ({/Times-Italic s})\"");
        fprintf(fp0, "set xrange [0:1000]\n");
        fprintf(fp0, "plot \"pg.dat\" with lines title \"P_g\" lc \"black\"\n");
        // subplot 2
        fprintf(fp0, "set size s1x*0.5, s1y*0.5\n");
        fprintf(fp0, "set origin o1x+s1x, o1y+s1y*0.1\n");
        fprintf(fp0, "%s\n", "set xlabel \"{/Times-Italic {v_p}} ({/Times-Italic km/s})\"");
        fprintf(fp0, "%s\n", "set ylabel \"{/Times-Italic depth}  ({/Times-Italic km})\"");
        fprintf(fp0, "set yrange [250:0]\n");
        fprintf(fp0, "set xrange [0:10]\n");
        fprintf(fp0, "set arrow from %lf,%lf  to %lf,%lf nohead lc \"red\"\n", vel[0].v, vel[0].d, vel[0].v, vel[1].d);
        fprintf(fp0, "set arrow from %lf,%lf  to %lf,%lf nohead lc \"red\"\n",
                vel[0].v, vel[1].d, vel[1].v, vel[1].d);
        fprintf(fp0, "set arrow from %lf,%lf  to %lf,%lf nohead lc \"red\"\n", vel[1].v, vel[1].d, vel[1].v, vel[2].d);
        fprintf(fp0, "set arrow from %lf,%lf  to %lf,%lf nohead lc \"red\"\n",
                vel[1].v, vel[2].d, vel[2].v, vel[2].d);
        fprintf(fp0, "set arrow from %lf,%lf  to %lf,%lf nohead lc \"red\"\n", vel[2].v, vel[2].d, vel[2].v, vel[3].d);
        fprintf(fp0, "set arrow from %lf,%lf  to %lf,%lf nohead lc \"red\"\n",
                vel[2].v, vel[3].d, vel[3].v, vel[3].d);
        fprintf(fp0, "set arrow from %lf,%lf  to %lf,%lf nohead lc \"red\"\n", vel[3].v, vel[3].d, vel[3].v, 250.);
        fprintf(fp0, "show arrow\n");
        fprintf(fp0, "set label \"hypocentre\" at 0.2,220\n");
        fprintf(fp0, "plot NaN t '', \"+\" using (8.5):(220) notitle lt 3 lc \"black\"\n");
        
        fclose(fp0);
    }
    else
        perror("pg2.sh");
    
    system("gnuplot pg.sh");
    
    system("gnuplot pg2.sh");
    system("open pg.pdf");
    system("open pg2.pdf");
    return 0;
}

int main()
{
    const int n = 5;
    Model * vel = (Model *)malloc(sizeof(Model) * n);
    vel[0].d = 0.0;
    vel[0].v = 4.0;
    vel[0].h = 15.0;
    
    vel[1].d = 15.0;
    vel[1].v = 6.0;
    vel[1].h = 50.0;
    
    vel[2].d = 65.0;
    vel[2].v = 8.0;
    vel[2].h = 100.0;
    
    vel[3].d = 165.0;
    vel[3].v = 8.5;
    vel[3].h = 200.0;
    
    vel[4].d = 365.0;
    vel[4].v = 8.75;
    vel[4].h = 500.0;
    
    
    double depth = 220;
    double v1 = 4;
    double v2 = 6;
    double h1 = 15;
    double x, t;
    double theta;
    double Delta = 1000;
    
    theta = get_take_off_angle(depth, v1, v2, h1, Delta, &t);
    printf("theta = %lf, traveltime = %lf, x = %lf\n", rad2deg(theta), t, x);
    
    theta = get_take_off_angle_Model(depth, vel, n, Delta, &t);
    printf("theta = %lf, traveltime = %lf, x = %lf\n", rad2deg(theta), t, x);
    plot_traveltime(vel, n, depth);
    free(vel);
    return 0;
}
