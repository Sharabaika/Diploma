#include "Solver.h"
#include <algorithm>
#include <iostream>

#define _USE_MATH_DEFINES

using std::find;
using std::min;
using std::max;

using std::cout;
using std::endl;

const double ONE_THIRD = 1.0 / 3.0;

const int CONDUCTOR_BORDER_INDEX = 1;
const int MEDIUM_REGION_INDEX = 2;
const int MEDIUM_OUTER_BORDER_INDEX = 3;

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

void SolveFast_Implementation(const Arr& x, const Arr& y, const JaggedArr& triangles, const ArrInt& segment_indices, const JaggedArr& trig_neighbors, const ArrInt& fluid_domain_nodes_indeces_array, const ArrInt& is_a_fluid_region_array, const ArrInt& is_a_wall_array, const double Pr, const double Ra, const double Ram, const double chi0, const double QPsi, const double QW, const double QT, const double max_error, const int max_cycles, Arr& Psi, Arr& W, Arr& T, const Arr& H_triangles, Arr& dHdx, Arr& dHdy, const Arr& normal_x, const Arr& normal_y, Arr& Delta_Psi, Arr& Delta_W, Arr& Delta_T)
{
    const int N_NODES = x.size();
    const int N_TRIANGLES = triangles.size();

    Delta_Psi = Arr();
    Delta_W = Arr();
    Delta_T = Arr();

    Arr Psi_new(N_NODES);
    Arr W_new(N_NODES);
    Arr T_new(N_NODES);

    int n_cycle = 0;
    double error = max_error * 2.0;

    Arr Psi_err(N_NODES);
    Arr W_err(N_NODES);
    Arr T_err(N_NODES);

    while (n_cycle < max_cycles && error >= max_error)
    {
        for (const int& n_node : fluid_domain_nodes_indeces_array)
        {
            double aPsi0 = 0.0;
            double aPsinb = 0.0;
            double aW0 = 0.0;
            double aWnb = 0.0;
            double aT0 = 0.0;
            double aTnb = 0.0;
            double S = 0.0;
            double I = 0.0;
            double source_integral = 0.0;

            const int& segment_index = segment_indices[n_node];

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

                if (!is_a_fluid_region_array[n0] || !is_a_fluid_region_array[n1] || !is_a_fluid_region_array[n2]) 
                {
                    continue;
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
                const double& T0 = T[n0]; const double& T1 = T[n1]; const double& T2 = T[n2];

                if (is_a_wall_array[n_node]) 
                {
                    // Boundaries # 
                    // ========== #

                    // Local Coords #
                    // ============ #

                    // mormal
                    const double sina = -normal_x[n_node];
                    const double cosa = normal_y[n_node];

                    const double& xc = x0; const double& yc = y0;
                    const double X0 = (x0 - xc) * cosa + (y0 - yc) * sina; const double Y0 = -(x0 - xc) * sina + (y0 - yc) * cosa;
                    const double X1 = (x1 - xc) * cosa + (y1 - yc) * sina; const double Y1 = -(x1 - xc) * sina + (y1 - yc) * cosa;
                    const double X2 = (x2 - xc) * cosa + (y2 - yc) * sina; const double Y2 = -(x2 - xc) * sina + (y2 - yc) * cosa;

                    double a = 0.0; double b = 0.0;
                    if (Y1 == 0.0) 
                    {
                        a = 0;
                        b = 2.0 * (Psi2 - Psi0) / Y2 * Y2;
                    }
                    else if (Y2 == 0) 
                    {
                        a = 0;
                        b = 2.0 * (Psi1 - Psi0) / Y1 * Y1;
                    }
                    else
                    {
                        const double DeltaGr = 0.5 * Y1 * Y2 * (X1 * Y2 - X2 * Y1);
                        const double DeltaAGr = 0.5 * ((Psi1 - Psi0) * Y2 * Y2 - (Psi2 - Psi0) * Y1 * Y1);
                        const double DeltaBGr = (Psi2 - Psi0) * X1 * Y1 - (Psi1 - Psi0) * X2 * Y2;

                        a = DeltaAGr / DeltaGr;
                        b = DeltaBGr / DeltaGr;
                    }

                    const double Xa = X1 * 0.5; const double Ya = Y1 * 0.5; 
                    const double Xb = X2 * 0.5; const double Yb = Y2 * 0.5;
                    const double Xc = (X1 + X2) * ONE_THIRD; const double Yc = (Y1 + Y2) * ONE_THIRD;

                    I = I + 0.5 * a * ((Yb * Yb - Ya * Ya) - (Xb * Xb - Xa * Xa));
                    I = I - 0.5 * b * ((Ya + Yc) * (Xc - Xa) + (Yb + Yc) * (Xb - Xc));

                    aWnb = aWnb + (7.0 * W1 + 7.0 * W2) / 216.0 * Delta;
                    aW0 = aW0 + 11.0 * Delta / 108.0;

                }
                else
                {
                    // Inner nodes #
                    // =========== #

                    // Stream function #
                    // =============== #
                    aPsi0 = aPsi0 + (x21 * x21 + y12 * y12) / Delta;
                    aPsinb = aPsinb + (Psi1 * (y20 * y12 + x02 * x21) + Psi2 * (y01 * y12 + x10 * x21)) / Delta;
                    S = S + (22.0 * W0 + 7.0 * W1 + 7.0 * W2) * Delta / 216.0;


                    // Stream local coords #
                    // =================== #
                    const double Apsi = Psi0 * y12 + Psi1 * y20 + Psi2 * y01;
                    const double Bpsi = Psi0 * x21 + Psi1 * x02 + Psi2 * x10;

                    const double Vx = Bpsi / Delta;
                    const double Vy = -Apsi / Delta;
                    const double Vel = sqrt(Vx * Vx + Vy * Vy);

                    double sina = 0.0;
                    double cosa = 1.0;

                    if (Vel > 1e-8) 
                    {
                        cosa = Vx / Vel;
                        sina = Vy / Vel;
                    }

                    const double xc = (x0 + x1 + x2) * ONE_THIRD;
                    const double yc = (y0 + y1 + y2) * ONE_THIRD;

                    const double X0 = (x0 - xc) * cosa + (y0 - yc) * sina;
                    const double Y0 = (y0 - yc) * cosa - (x0 - xc) * sina;
                    const double X1 = (x1 - xc) * cosa + (y1 - yc) * sina;
                    const double Y1 = (y1 - yc) * cosa - (x1 - xc) * sina;
                    const double X2 = (x2 - xc) * cosa + (y2 - yc) * sina;
                    const double Y2 = (y2 - yc) * cosa - (x2 - xc) * sina;

                    const double Xmax = max({ X0, X1, X2 });
                    const double Xmin = min({ X0, X1, X2 });

                    const double X12 = X1 - X2; const double Y12 = Y1 - Y2;
                    const double Y20 = Y2 - Y0; const double Y01 = Y0 - Y1;
                    const double X10 = X1 - X0; const double X02 = X0 - X2;

                    // Temperature //
                    // =========== //
                    const double abT = Vel * Y12 * Y0 / 8.0 + X12 / 2.0;
                    const double acT = Vel * Y12 / 2.0;

                    double k0T = 0;
                    double k1T = 0;
                    double k2T = 0;

                    if (abs(Vel) < 1e-14)
                    {

                        k0T = 0;
                        k1T = 0;
                        k2T = 0;
                    }
                    else if (abs(Vel * (Xmax - Xmin)) < 1e-8)
                    {
                        const double DT = Vel * (X0 * Y12 + X1 * Y20 + X2 * Y01);
                        k0T = (abT * Vel * (X2 - X1) + acT * (-Y12 + Vel * ((X1 - Xmax) * Y2 - (X2 - Xmax) * Y1))) / DT;
                        k1T = (abT * Vel * (X0 - X2) + acT * (-Y20 + Vel * ((X2 - Xmax) * Y0 - (X0 - Xmax) * Y2))) / DT;
                        k2T = (abT * Vel * (X1 - X0) + acT * (-Y01 + Vel * ((X0 - Xmax) * Y1 - (X1 - Xmax) * Y0))) / DT;
                    }
                    else
                    {
                        const double E0T = exp(Vel * (X0 - Xmax));
                        const double E1T = exp(Vel * (X1 - Xmax));
                        const double E2T = exp(Vel * (X2 - Xmax));
                        const double DT = E0T * Y12 + E1T * Y20 + E2T * Y01;
                        k0T = (abT * (E2T - E1T) + acT * (E1T * Y2 - E2T * Y1)) / DT;
                        k1T = (abT * (E0T - E2T) + acT * (E2T * Y0 - E0T * Y2)) / DT;
                        k2T = (abT * (E1T - E0T) + acT * (E0T * Y1 - E1T * Y0)) / DT;
                    }

                    const double AT = (T0 * y12 + T1 * y20 + T2 * y01) / Delta;
                    const double BT = (T0 * x21 + T1 * x02 + T2 * x10) / Delta;

                    const double c_mag = dHdx[n_trig_neighbour] * BT - dHdy[n_trig_neighbour] * AT;
                    const double mu = chi0 * H_triangles[n_trig_neighbour] / (1 + chi0 * H_triangles[n_trig_neighbour]);
                    source_integral += (Ra * AT + Ram * mu * c_mag) * Delta / 6.0;

                    aT0 = aT0 + k0T;
                    aTnb = aTnb + k1T * T1 + k2T * T2;


                    // Vorticity //
                    // ========= //

                    const double VelPr = Vel / Pr;
                    const double abW = VelPr * Y12 * Y0 / 8.0 + X12 / 2.0;
                    const double acW = VelPr * Y12 / 2.0;

                    double k0W = 0;
                    double k1W = 0;
                    double k2W = 0;

                    if (abs(Vel) < 1e-14)
                    {
                        k0W = 0;
                        k1W = 0;
                        k2W = 0;
                    }
                    else if (abs(VelPr * (Xmax - Xmin)) < 1e-8)
                    {
                        const double DW = VelPr * (X0 * Y12 + X1 * Y20 + X2 * Y01);
                        k0W = (abW * VelPr * (X2 - X1) + acW * (-Y12 + VelPr * ((X1 - Xmax) * Y2 - (X2 - Xmax) * Y1))) / DW;
                        k1W = (abW * VelPr * (X0 - X2) + acW * (-Y20 + VelPr * ((X2 - Xmax) * Y0 - (X0 - Xmax) * Y2))) / DW;
                        k2W = (abW * VelPr * (X1 - X0) + acW * (-Y01 + VelPr * ((X0 - Xmax) * Y1 - (X1 - Xmax) * Y0))) / DW;
                    }
                    else
                    {
                        const double E0 = exp(VelPr * (X0 - Xmax));
                        const double E1 = exp(VelPr * (X1 - Xmax));
                        const double E2 = exp(VelPr * (X2 - Xmax));
                        const double DW = (E0 * Y12 + E1 * Y20 + E2 * Y01);
                        k0W = (abW * (E2 - E1) + acW * (E1 * Y2 - E2 * Y1)) / DW;
                        k1W = (abW * (E0 - E2) + acW * (E2 * Y0 - E0 * Y2)) / DW;
                        k2W = (abW * (E1 - E0) + acW * (E0 * Y1 - E1 * Y0)) / DW;
                    }

                    aW0 = aW0 + k0W;
                    aWnb = aWnb + k1W * W1 + k2W * W2;
                }

                // Advance //
                // ======= //

                if (!is_a_wall_array[n_node]) 
                {
                    // Medium //
                    // ====== //

                    Psi_new[n_node] = (-aPsinb + S) / aPsi0;
                    W_new[n_node] = -(aWnb + source_integral) / aW0;
                    T_new[n_node] = -aTnb / aT0;
                }
                else
                {
                    // Wall //
                    // ==== //

                    Psi_new[n_node] = 0;
                    W_new[n_node] = -(I + aWnb) / aW0;
                    if (segment_index == CONDUCTOR_BORDER_INDEX)
                    {
                        T_new[n_node] = 1;
                    }
                    else
                    {
                        T_new[n_node] = 0;
                    }
                }

                Psi_err[n_node] = (Psi[n_node] - Psi_new[n_node]);
                W_err[n_node] = (W_new[n_node] - W[n_node]);
                T_err[n_node] = (T_new[n_node] - T[n_node]);
            }
        }

        double Psi_sum = 0, W_sum = 0, T_sum = 0;
        double Psi_errors = 0, W_errors = 0, T_errors = 0;
        for (int n_node = 0; n_node < N_NODES; n_node++)
        {
            Psi_sum += Psi_new[n_node] * Psi_new[n_node];
            Psi_errors += Psi_err[n_node] * Psi_err[n_node];
            Psi[n_node] = Psi[n_node] * (1.0 - QPsi) + Psi_new[n_node] * QPsi;

            W_sum += W_new[n_node] * W_new[n_node];
            W_errors += W_err[n_node] * W_err[n_node];
            W[n_node] = W[n_node] * (1.0 - QW) + W_new[n_node] * QW;

            T_sum += T_new[n_node] * T_new[n_node];
            T_errors += T_err[n_node] * T_err[n_node];
            T[n_node] = T[n_node] * (1.0 - QT) + T_new[n_node] * QT;
        }

        const double Delta_Psi_Error = sqrt(Psi_errors / Psi_sum) / QPsi;
        const double Delta_W_Error = sqrt(W_errors / W_sum) / QW;
        const double Delta_T_Error = sqrt(T_errors / T_sum) / QT;

        if (n_cycle % 5000 == 0)
        {
            printf("n cycle = %06i, dw = %.5e, dpsi = %.5e, dt = %.5e\n", n_cycle, Delta_W_Error, Delta_Psi_Error, Delta_T_Error);
        }

        Delta_Psi.push_back(Delta_Psi_Error);
        Delta_W.push_back(Delta_W_Error);
        Delta_T.push_back(Delta_T_Error);

        error = max({ Delta_Psi_Error, Delta_W_Error, Delta_T_Error});

        n_cycle++;
    }
}

void SolveMagnetics_fast(const Arr& x, const Arr& y, const JaggedArr& triangles, const ArrInt& segment_indices, const JaggedArr& trig_neighbors, const ArrInt& trianlge_indices, const double chi0, const double H0, const double mu0, const double QFi, const double max_error, const int max_cycles, Arr& H, Arr& Mu, Arr& Fi, Arr& H_nodes, Arr& Delta_Fi)
{
    const int N_NODES = x.size();
    const int N_TRIANGLES = triangles.size();

    Delta_Fi = Arr();

    Arr Fi_new(N_NODES);
    Arr H_new(N_TRIANGLES);
    Arr Mu_new(N_TRIANGLES);
    
    H_nodes = Arr(N_NODES);

    Arr Fi_err(N_NODES);
    
    double error = 2.0 * max_error;
    int n_cycle = 0;
    while (100 > n_cycle || (n_cycle < max_cycles && error >= max_error))
    {
        for (int n_node = 0; n_node < N_NODES; n_node++)
        {
            double a0F = 0.0;
            double anbF = 0.0;


            for (int n_trig_neighbour : trig_neighbors[n_node])
            {
                const ArrInt Triangle = triangles[n_trig_neighbour];

                int n0 = Triangle[0];
                int n1 = Triangle[1];
                int n2 = Triangle[2];
                if (n_node == n1)
                {
                    n0 = n_node;
                    n1 = Triangle[2];
                    n2 = Triangle[0];
                }
                else if (n_node == n2)
                {
                    n0 = n_node;
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


                // Fi //
                // == //

                a0F = a0F + Mu[n_trig_neighbour] * 0.5 / Delta * (y12 * y12 + x21 * x21);
                anbF = anbF + Mu[n_trig_neighbour] * 0.5 / Delta * (Fi[n1] * (y12 * y20 + x21 * x02) + Fi[n2] * (y12 * y01 + x21 * x10));
            }

            const int segment_index = segment_indices[n_node];
            const bool bInfinity = segment_index == 5;
            if (!bInfinity) {
                if (abs(a0F) > 0)
                {
                    Fi_new[n_node] = -anbF / a0F;
                }
            }
            else
            {
                Fi_new[n_node] = H0 * y[n_node];
            }

            Fi_err[n_node] = Fi_new[n_node] - Fi[n_node];
        }

        for (int n_triangle = 0; n_triangle < N_TRIANGLES; n_triangle++)
        {
            const auto& Triangle = triangles[n_triangle];
            int n0 = Triangle[0];
            int n1 = Triangle[1];
            int n2 = Triangle[2];

            const double& x0 = x[n0]; const double& y0 = y[n0];
            const double& x1 = x[n1]; const double& y1 = y[n1];
            const double& x2 = x[n2]; const double& y2 = y[n2];

            const double x10 = x1 - x0; const double y01 = y0 - y1;
            const double x21 = x2 - x1; const double y12 = y1 - y2;
            const double x02 = x0 - x2; const double y20 = y2 - y0;
            const double Delta = x10 * y20 - x02 * y01;


            const double DeltaA = Fi[n0] * y12 + Fi[n1] * y20 + Fi[n2] * y01;
            const double DeltaB = Fi[n0] * x21 + Fi[n1] * x02 + Fi[n2] * x10;

            // H //
            // = //
            const double Hx = DeltaA / Delta;
            const double Hy = DeltaB / Delta;
            H_new[n_triangle] = sqrt(Hx * Hx + Hy * Hy);


            // Mu //
            // == //
            const int segment_index = trianlge_indices[n_triangle];
            if (segment_index == 0)
            {
                Mu_new[n_triangle] = mu0;
            }
            else if (segment_index == 2)
            {
                Mu_new[n_triangle] = 1.0 + chi0 / (1.0 + chi0 * H_new[n_triangle]);
            }
            else if (segment_index == 4)
            {
                Mu_new[n_triangle] = 1.0;
            }
        }

        Mu = Mu_new;
        H = H_new;

        double numerator = 0.0;
        double denomenator = 0.0;
        for (int n_node = 0; n_node < N_NODES; n_node++)
        {
            numerator += Fi_err[n_node] * Fi_err[n_node];
            denomenator += Fi[n_node] * Fi[n_node];

            Fi[n_node] = Fi[n_node] * (1.0 - QFi) + Fi_new[n_node] * QFi;
        }
        error = sqrt(numerator / denomenator) / QFi;

        if (n_cycle % 50 == 0)
        {
            printf("n cycle = %06i, dFi = %.5e \n", n_cycle, error);
        }

        Delta_Fi.push_back(error);
        n_cycle++;
    }


    for (size_t n_node = 0; n_node < N_NODES; n_node++)
    {
        double numerator = 0.0;
        double denomenator = 0.0;
        for (int n_trig_neighbour : trig_neighbors[n_node])
        {
            const ArrInt Triangle = triangles[n_trig_neighbour];

            int n0 = Triangle[0];
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

            const double Mx = (x1 + x2) * 0.5;
            const double My = (y1 + y2) * 0.5;

            const double dx = Mx - x0;
            const double dy = My - y0;
            const double mr = sqrt(dx * dx + dy * dy);
            const double r = 2.0 * mr / 3.0;

            numerator += H[n_trig_neighbour] / r;
            denomenator += 1.0 / r;
        }
        
        H_nodes[n_node] = numerator / denomenator;
       
    }

}
