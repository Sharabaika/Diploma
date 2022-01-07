#include "Solver.h"
#include <algorithm>
#include <iostream>

#define _USE_MATH_DEFINES

using std::find;
using std::min;
using std::max;

using std::cout;
using std::endl;

bool containts_index(const ArrInt& indecies, const int& index) 
{
    return find(indecies.begin(), indecies.end(), index) != indecies.end();
}

void SolveFluid_Implementation(const Arr x,const Arr y,const JaggedArr triangles,const JaggedArr segment_indices,const JaggedArr trig_neighbors,const JaggedArr node_neighbours,
	const double Re,const double Vx, 
	const double QPsi,const double QW,const double max_error,const int max_cycles,
	Arr& Psi, Arr& W,
	Arr& Delta_Psi, Arr& Delta_W
)
{
    const clock_t begin_time = std::clock();

    int n_cycle = 0;
    int MAX_CYCLES = max_cycles;

    double error = max_error * 2.0;

    const int N_NODES = x.size();

    const double Pr = 1.0 / Re;
    const double Vc = 1.0;

    Psi = Arr(N_NODES);
    W = Arr(N_NODES);

    Delta_Psi = Arr();
    Delta_W = Arr();

    Arr Psi_new(N_NODES);
    Arr W_new(N_NODES);


    for (int n_node = 0; n_node < N_NODES; n_node++)
    {
        const ArrInt& Segment_Index = segment_indices[n_node];
        const bool b_wall = containts_index(Segment_Index, SOLID_BORDER_INDEX);
        if (b_wall) 
        {
            Psi[n_node] = 0;
        }
        else
        {
            const double& xn = x[n_node];
            const double& yn = y[n_node];

            Psi[n_node] = 0.0001 * sin(xn * PI) * sin(yn * PI);
        }

        W[n_node] = 0.0;
    }

    while (n_cycle < max_cycles && error >= max_error)
    { 
        for (int n_node = 0; n_node < N_NODES; n_node++)
        {
            const ArrInt& Segment_Index = segment_indices[n_node];

            // USUAL MEDIUM //
            // ============ //
            if (containts_index(Segment_Index, MEDIUM_INDEX))
            {
                double Psi_BorderIntegral_a0 = 0;
                double Psi_BorderIntegral_nb = 0;
                double Psi_AreaIntegral = 0;

                double W_BorderIntegral = 0;
                double W_BorderIntegral_k0 = 0;
                double W_AreaIntegral = 0;

                for (const int& n_trig_neighbour : trig_neighbors[n_node])
                {
                    const ArrInt Triangle = triangles[n_trig_neighbour];

                    int n0 = n_node;
                    int n1 = Triangle[1];
                    int n2 = Triangle[2];
                    if (n_node == n1)
                    {
                        n1 = Triangle[2];
                        n2 = Triangle[0];
                    }
                    else if (n_node == n2)
                    {
                        n1 = Triangle[0];
                        n2 = Triangle[1];
                    }

                    const double& x0 = x[n0]; const double& y0 = y[n0];
                    const double& x1 = x[n1]; const double& y1 = y[n1];
                    const double& x2 = x[n2]; const double& y2 = y[n2];

                    const double x10 = x1 - x0; const double y01 = y0 - y1;
                    const double x21 = x2 - x1; const double y12 = y1 - y2;
                    const double x02 = x0 - x2; const double y20 = y2 - y0;
                    const double Delta = x10 * y20 - x02 * y01;

                    const double& Psi0 = Psi[n0]; const double& Psi1 = Psi[n1]; const double& Psi2 = Psi[n2];
                    const double& W0 = W[n0]; const double& W1 = W[n1]; const double& W2 = W[n2];


                    // Psi //
                    // === //
                    const double Delta_PsiA = Psi0 * y12 + Psi1 * y20 + Psi2 * y01;
                    const double Delta_PsiB = Psi0 * x21 + Psi1 * x02 + Psi2 * x10;

                    const double A_Psi = Delta_PsiA / Delta;
                    const double B_PSi = Delta_PsiB / Delta;

                    const double a0 = (x21 * x21 + y12 * y12) / Delta;
                    const double a1_psi1 = Psi1 * (y20 * y12 + x02 * x21) / Delta;
                    const double a2_psi2 = Psi2 * (y01 * y12 + x10 * x21) / Delta;

                    Psi_BorderIntegral_a0 += a0;
                    Psi_BorderIntegral_nb += (a1_psi1 + a2_psi2);

                    Psi_AreaIntegral += (22.0 * W0 + 7.0 * W1 + 7.0 * W2) * Delta / 216.0;

                    // W //
                    // = //
                    const double U_x = B_PSi;
                    const double U_y = -A_Psi;
                    const double U_ell = sqrt(U_x * U_x + U_y * U_y);

                    double sina = 0;
                    double cosa = 1;

                    if (U_ell > 0.0)
                    {
                        sina = U_y / U_ell;
                        cosa = U_x / U_ell;
                    }

                    const double xc = (x0 + x1 + x2) / 3.0;
                    const double yc = (y0 + y1 + y2) / 3.0;

                    const double X0 = (x0 - xc) * cosa + (y0 - yc) * sina; const double Y0 = -(x0 - xc) * sina + (y0 - yc) * cosa;
                    const double X1 = (x1 - xc) * cosa + (y1 - yc) * sina; const double Y1 = -(x1 - xc) * sina + (y1 - yc) * cosa;
                    const double X2 = (x2 - xc) * cosa + (y2 - yc) * sina; const double Y2 = -(x2 - xc) * sina + (y2 - yc) * cosa;

                    const double X12 = X1 - X2; const double Y12 = Y1 - Y2;
                    const double X20 = X2 - X0; const double Y20 = Y2 - Y0;
                    const double X01 = X0 - X1; const double Y01 = Y0 - Y1;

                    const double X_min = min({ X0, X1, X2 });
                    const double X_max = max({ X0, X1, X2 });

                    double kw0 = 0; double kw1 = 0; double kw2 = 0;

                    const double ReVel = U_ell / (Pr * Vc);
                    const double vel = ReVel * (X_max - X_min);

                    const double aBw = U_ell * Y12 * Y0 / (8.0 * Pr) + Vc * X12 / 2.0;
                    const double aCw = (U_ell * Y12) / (2.0 * Pr);

                    if (vel <= 1e-8)
                    {
                        const double DW = ReVel * (X0 * Y12 + X1 * Y20 + X2 * Y01);

                        kw0 = (aBw * ReVel * (X2 - X1) + aCw * (-Y12 + ReVel * ((X1 - X_max) * Y2 - (X2 - X_max) * Y1))) / DW;
                        kw1 = (aBw * ReVel * (X0 - X2) + aCw * (-Y20 + ReVel * ((X2 - X_max) * Y0 - (X0 - X_max) * Y2))) / DW;
                        kw2 = (aBw * ReVel * (X1 - X0) + aCw * (-Y01 + ReVel * ((X0 - X_max) * Y1 - (X1 - X_max) * Y0))) / DW;
                    }
                    else
                    {
                        const double EW0 = exp(U_ell * (X0 - X_max) / (Pr * Vc));
                        const double EW1 = exp(U_ell * (X1 - X_max) / (Pr * Vc));
                        const double EW2 = exp(U_ell * (X2 - X_max) / (Pr * Vc));

                        const double Delta_W = EW0 * Y12 + EW1 * Y20 + EW2 * Y01;
                        const double one_over_Delta_W = 1.0 / Delta_W;

                        kw0 = (aBw * (EW2 - EW1) + aCw * (EW1 * Y2 - EW2 * Y1)) * one_over_Delta_W;
                        kw1 = (aBw * (EW0 - EW2) + aCw * (EW2 * Y0 - EW0 * Y2)) * one_over_Delta_W;
                        kw2 = (aBw * (EW1 - EW0) + aCw * (EW0 * Y1 - EW1 * Y0)) * one_over_Delta_W;
                    }

                    W_BorderIntegral += W1 * kw1 + W2 * kw2;
                    W_BorderIntegral_k0 += kw0;
                }

                Psi_new[n_node] = (-Psi_BorderIntegral_nb + Psi_AreaIntegral) / Psi_BorderIntegral_a0;
                W_new[n_node] = (-W_BorderIntegral + W_AreaIntegral) / W_BorderIntegral_k0;

            }
            // SOLID BORDER //
            // ============ //
            else if (containts_index(Segment_Index, SOLID_BORDER_INDEX)) 
            {
                // normal //
                // ------ //
                 
                double sina = 0, cosa = 0;

                if (containts_index(Segment_Index, 10)) 
                {
                    // left
                    sina = -1; cosa = 0;
                }
                else if (containts_index(Segment_Index, 11)) 
                {
                    // right
                    sina = 1; cosa = 0;
                }
                else if (containts_index(Segment_Index, 12)) 
                {
                    // bottom
                    sina = 0; cosa = 1;
                }
                else if (containts_index(Segment_Index, 13)) 
                {
                    // upp
                    sina = 0; cosa = -1;
                }

                // local coords // 
                // ------------ //
                const double xc = x[n_node];
                const double yc = y[n_node];

                // Integral values //
                // --------------- //
                // Just an area of triangle, nothing more
                double W_Source_Integral = 0;

                // Integral of W over triangle
                double W_Source_Area_Integral = 0;

                // Integral a to b
                double W_Border_Integral = 0;

                for (const int& n_trig_neighbour : trig_neighbors[n_node])
                {
                    const ArrInt Triangle = triangles[n_trig_neighbour];

                    int n0 = n_node;
                    int n1 = Triangle[1];
                    int n2 = Triangle[2];
                    if (n_node == n1)
                    {
                        n1 = Triangle[2];
                        n2 = Triangle[0];
                    }
                    else if (n_node == n2)
                    {
                        n1 = Triangle[0];
                        n2 = Triangle[1];
                    }

                    const double& x0 = x[n0]; const double& y0 = y[n0];
                    const double& x1 = x[n1]; const double& y1 = y[n1];
                    const double& x2 = x[n2]; const double& y2 = y[n2];

                    const double& Psi0 = Psi[n0]; const double& Psi1 = Psi[n1]; const double& Psi2 = Psi[n2];
                    const double& W0 = W[n0]; const double& W1 = W[n1]; const double& W2 = W[n2];

                    const double X0 = (x0 - xc) * cosa + (y0 - yc) * sina; const double Y0 = -(x0 - xc) * sina + (y0 - yc) * cosa;
                    const double X1 = (x1 - xc) * cosa + (y1 - yc) * sina; const double Y1 = -(x1 - xc) * sina + (y1 - yc) * cosa;
                    const double X2 = (x2 - xc) * cosa + (y2 - yc) * sina; const double Y2 = -(x2 - xc) * sina + (y2 - yc) * cosa;

                    bool b_moving_border = containts_index(Segment_Index, 13);
                    const double v_border_local = b_moving_border ? Vx : 0.0;

                    const double c1 = Psi1 - Psi0 - v_border_local * Y1;
                    const double c2 = Psi2 - Psi0 - v_border_local * Y2;

                    const double delta = 0.5 * Y1 * Y2 * (X1 * Y2 - X2 * Y1);
                    const double delta_a = 0.5 * (c1 * (Y2 * Y2) - c2 * (Y1 * Y1));
                    const double delta_b = c2 * X1 * Y1 - c1 * X2 * Y2;

                    double a = 0, b = 0;
                    if (delta != 0)
                    {
                        a = delta_a / delta;
                        b = delta_b / delta;
                    }
                    else
                    {
                        if (Y1 == 0)
                        {
                            b = 2 * c2 / (Y2 * Y2);
                        }
                        else if (Y2 == 0)
                        {
                            b = 2 * c1 / (Y1 * Y1);
                        }
                    }

                    const double xa = X1 * 0.5; const double ya = Y1 * 0.5;
                    const double xb = X2 * 0.5; const double yb = Y2 * 0.5;
                    const double xc = (X1 + X2) / 3.0; const double yc = (Y1 + Y2) / 3.0;

                    const double moving_border_add = v_border_local * (xb - xa);

                    const double border_part = 0.5 * a * (yb * yb - ya * ya) - 0.5 * a * (xb * xb - xa * xa) - 
                        b * (0.5 * (ya + yc) * (xc - xa) + 0.5 * (yb + yc) * (xb - xc)) + 
                        moving_border_add;

                    const double x10 = x1 - x0; const double y01 = y0 - y1;
                    const double x21 = x2 - x1; const double y12 = y1 - y2;
                    const double x02 = x0 - x2; const double y20 = y2 - y0;
                    const double triangle_delta = x10 * y20 - x02 * y01;

                    W_Border_Integral += border_part;
                    W_Source_Area_Integral += 11.0 * triangle_delta / 108.0;
                    W_Source_Integral += (7.0 * W1 + 7.0 * W2) * triangle_delta / 216.0;
                }

                Psi_new[n_node] = 0.0;

                double xn = x[n_node];
                double yn = y[n_node];

                if (xn == 0 && yn == 0 ||
                    xn == 0 && yn == 1 ||
                    xn == 1 && yn == 1 ||
                    xn == 1 && yn == 0) 
                {
                    W_new[n_node] = 0;
                }
                else
                {
                    W_new[n_node] = (-W_Border_Integral + W_Source_Integral) / W_Source_Area_Integral;
                }
            }
        }



        double Psi_sum = 0, W_sum = 0;
        double Psi_errors = 0, W_errors = 0;
        for (int n_node = 0; n_node < N_NODES; n_node++)
        {
            const double& psi_old = Psi[n_node];
            const double& psi_new = Psi_new[n_node];

            if (psi_new != 0)
            {
                Psi_sum += psi_new * psi_new;
                Psi_errors += (psi_old - psi_new) * (psi_old - psi_new);
            }

            Psi[n_node] = psi_new * QPsi + psi_old * (1.0 - QPsi);

            const double& w_old = W[n_node];
            const double& w_new = W_new[n_node];

            if (w_new != 0)
            {
                W_sum = w_new * w_new;
                W_errors += (w_old - w_new) * (w_old - w_new);
            }

            W[n_node] = w_new * QW + w_old * (1.0 - QW);
        }

        const double Delta_Psi_Error = sqrt(Psi_errors/Psi_sum) / QPsi;
        const double Delta_W_Error = sqrt(W_errors / W_sum) / QW;

        if (n_cycle % 50 == 0) 
        {
            printf("n cycle = %06i, dw = %.5e, dpsi = %.5e \n", n_cycle, Delta_W_Error, Delta_Psi_Error);
        }

        Delta_Psi.push_back(Delta_Psi_Error);
        Delta_W.push_back(Delta_W_Error);

        error = max({ Delta_Psi_Error, Delta_W_Error });

        n_cycle++;
    }

    const clock_t end_time = clock();

    printf("===================================================\n");
    printf("n cycle = %06i, dw = %.5e, dpsi = %.5e \n", n_cycle, Delta_W[Delta_W.size() - 1], Delta_Psi[Delta_Psi.size() - 1]);
    printf("time spent = %f \n", float(end_time - begin_time) / CLOCKS_PER_SEC);
}