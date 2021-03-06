
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <getopt.h>

int factorial(int n);

int main(int argc, char **argv){

    // dynamic variables
    double xmin   =   -2.5;            // minimal x-value, set with -a or --xmin
    double xmax   =    2.5;            // maximal x-value, set with -b or --xmax
    double dx     =    0.025;          // stepsize of x,   set with -d or --x-spacing
    double k      = 4594.109225112563; // force-constant in kJ/mol/angstrom^2, set with -k or --force-constant
    double mu     =    1.0;            // reduced mass of involved particles in g/mol, set with -m or --reduced-mass
    char *outputfile = "/dev/stdout";
    static int kcal_flag = 0;          // set flag to kcal/mol
    static int evecs_flag = 0;
    int numberofeigenstates=6;

    int c;
    while(1){
        static struct option long_options[] = {
            /* These options set a flag. */
            {"kcal" , no_argument, &kcal_flag , 1},
            {"evecs", no_argument, &evecs_flag, 1},
            /* These options don’t set a flag.
               We distinguish them by their indices. */
            {"help",                  no_argument, 0, 'h'},
            {"xmin",                  required_argument, 0, 'a'},
            {"xmax",                  required_argument, 0, 'b'},
            {"x-spacing",             required_argument, 0, 'd'},
            {"force-constant",        required_argument, 0, 'k'},
            {"reduced-mass",          required_argument, 0, 'm'},
            {"output-file",           required_argument, 0, 'o'},
            {"number-of-eigenstates", required_argument, 0, 'n'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "a:b:d:k:m:o:n:h", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch(c){
            case 0:
            /* If this option set a flag, do nothing else now. */
            if (long_options[option_index].flag != 0)
                break;
            printf("option %s", long_options[option_index].name);
            if (optarg)
                printf(" with arg %s", optarg);
            printf("\n");
            break;

            case 'h':
                printf("%s\t", argv[0]);
                printf("[--help|-h]");
                printf("[--xmin|-a <value>]");
                printf("[--xmax|-b <value>]");
                printf("\n\t\t\t");
                printf("[--x-spacing|-d <value>]");
                printf("[--force-constant|-k <value>]");
                printf("\n\t\t\t");
                printf("[--reduced-mass|-m <value>]");
                printf("[--output-file|-o]");
                printf("[--kcal]");
                printf("\n\t\t\t");
                printf("[--number-of-eigenstates|-n]");
                printf("\n\n");

                printf("-h, --help\t\t\tShow this help dialogue\n");
                printf("-a, --xmin\t\t\tSet minimum x-value\n");
                printf("-b, --xmax\t\t\tSet maximum x-value\n");
                printf("-d, --x-spacing\t\t\tSet spacing between x-values\n");
                printf("-k, --force-constant\t\tSet force constant in kJ/(mol*angstrom^2)\n");
                printf("-m, --reduced-mass\t\tSet reduced mass of involved particles in g/mol\n");
                printf("-o, --output-file\t\tName of output file\n");
                printf("-n, --number-of-eigenstates\tNumber of calculated eigenstates\n");
                printf("    --kcal\t\t\tOutput in kcal/mol\n");
                printf("    --evecs\t\t\tOutput eigenvectors directly (don't add evals)");
                printf("\n");
                exit (0);

            case 'a':
                //printf ("option -a with value `%s'\n", optarg);
                xmin = atof(optarg);
                break;

            case 'b':
                xmax = atof(optarg);
                break;

            case 'd':
                dx = atof(optarg);
                break;

            case 'k':
                k = atof(optarg);
                break;

            case 'm':
                mu = atof(optarg);
                break;

            case 'o':
                outputfile = optarg;
                break;

            case 'n':
                numberofeigenstates = atoi(optarg);
                break;

            default:
                abort();
        }
    }

    int i;
    double x;
    double omega, term1, term2, term3, arg, eval;

    double H[numberofeigenstates];
    double psi[numberofeigenstates];
    FILE *fd = fopen(outputfile, "w");

    // constants
    double pi       = 3.1415926535897932384626433832795;
    double planck   = 6.626070040E-34;      // Js
    double avogadro = 6.022140857E23;       // 1/mol
    double amu      = 1/(1000*avogadro);    // kg/(g/mol) <=> 1/1000/avogadro

    // Psi(n,x) = (mu*omega/pi/hbar)^(1/4)      // term1
    //            * 1/sqrt(n!*2^n)*
    //            * H_n(sqrt(mu*omega/hbar)*x)  // H_n(term2*x)
    //            * exp(-mu*omega*x^2/2/hbar)   // exp(term3*x^2)

    // omega^2 = k kJ/mol/angstrom^2 * 1/mu mol/g * 1000 kg*m^2/s^2/kJ * 1000 g/kg * 10^20 angstrom^2/m^2
    //         = 10^26 * k/mu 1/s^2
    omega = 1.0E13*sqrt(k/mu);  // 1/s^2

    // term1^4 = mu g/mol * omega 1/s * hbar 1/(J*s) * amu kg/(g/mol) / 10^20 m^2/angstrom^2
    //         = 10^-20*amu * mu*omega/pi/hbar 1/angstrom^2
    term1 = 1.0E-5*pow(2.0*mu*amu*omega/planck, 0.25); // 1/sqrt(angstrom)

    // term2^2 = mu g/mol * omega 1/s * hbar 1/(J*s) * amu kg/(g/mol) / 10^20 m^2/angstrom^2
    //         = 10^-20*amu * mu*omega/hbar 1/angstrom^2
    term2 = 1.0E-10*sqrt(mu*amu*omega*2.0*pi/planck); // 1/angstrom

    // term3   = mu g/mol * omega 1/s * hbar 1/(J*s) * amu kg/(g/mol) / 10^20 m^2/angstrom^2
    //         = 10^-20*amu * -mu*omega/2/hbar 1/angstrom^2
    term3 = -1.0E-20*mu*amu*omega*pi/planck;         // 1/angstrom^2

    // analytical eigenvalues
    // eval    = hbar J*s * omega 1/s * avogadro 1/mol / 1000 kJ/J
    eval = planck/pi*omega*avogadro/2000.0;         // kJ/mol

    if(kcal_flag == 1){
        eval = eval/4.184;
    }

    fprintf(fd, "# Harmonic Oscillator\n");
    fprintf(fd, "#  V(x) = 0.5*k*(x-x0)^2 = 0.5*mu*omega^2*(x-x0)^2\n");
    fprintf(fd, "# omega:            %+18.12e 1/s\n", omega);
    fprintf(fd, "# reduced mass:     %+18.12e g/mol\n", mu);
    fprintf(fd, "# prefactor:        %+18.12e 1/sqrt(angstrom)\n", term1);
    fprintf(fd, "# hermite-argument: %+18.12e 1/angstrom\n", term2);
    fprintf(fd, "# exponent:         %+18.12e 1/angstrom^2\n", term3);

    if(kcal_flag == 1)  fprintf(fd, "# Energies in kcal/mol:\n");
    else                fprintf(fd, "# Energies in kJ/mol:\n");

    for(i=0; i<numberofeigenstates; ++i){
            fprintf(fd, "#\tE%02d = % 26.20lf\n", i, (0.5+i)*eval);
    }

    // Calculate E+Psi and output data
    fprintf(fd, "#              x                        V             ");
    for(i=0; i<numberofeigenstates; ++i){
        fprintf(fd, "         E%02d + Psi%02d        ", i, i);
    }
    fprintf(fd, "\n");
    for(x = xmin; x <= xmax; x += dx){

        arg = term2*x;

        // H(n,x) = (-1)^n * e^(x^2) * d^n/dx^n e^(-x^2)
        H[0] = 1;
        H[1] = 2*arg;
        for(i=2; i<numberofeigenstates; ++i){
            H[i] = 2.0*arg*H[i-1] - 2.0*((double)i - 1.0)*H[i-2];
        }
        for(i=0; i<numberofeigenstates; ++i){
            psi[i] = term1 * 1.0/sqrt(pow(2,i) * (double)factorial(i)) * H[i] * exp(term3*x*x);
        }

        fprintf(fd,  "% 26.20lf  ", x);
        if(kcal_flag == 1){
            fprintf(fd, "% 26.20lf  ", 0.5*k*x*x/4.184);
        }else{
            fprintf(fd, "% 26.20lf  ", 0.5*k*x*x);
        }

        if(evecs_flag == 0){
            for(i=0; i<numberofeigenstates; ++i){
                fprintf(fd, "% 26.20lf  ", psi[i] + (0.5+(double)i)*eval);
            }
            fprintf(fd, "\n");
        } else{
            for(i=0; i<numberofeigenstates; ++i){
                fprintf(fd, "% 26.20lf  ", psi[i]);
            }
            fprintf(fd, "\n");
        }
    }
    fclose(fd);

    return 0;
}

int factorial(int n){
    int fact = 1;

    while(n > 1){
        fact = fact * n;
        n = n - 1;
    }
    return fact;
}
