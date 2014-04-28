#ifndef POLYWOG_LOGIT_H
#define POLYWOG_LOGIT_H

inline
double transformLogit(double x)
{
    return 1 / (1 + exp(-1.0 * x));
}

#endif  // POLYWOG_LOGIT_H
