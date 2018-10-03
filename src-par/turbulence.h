#ifndef TURB_H
#define TURB_H

real ***filter_( real ***U, real h[3]);
void DynamicSmagorinsky(real ***rho, real ***u, real ***v, real ***w, real ***mu_SGS);

#endif