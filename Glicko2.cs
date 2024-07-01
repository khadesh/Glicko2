using System;

namespace EloCalculations
{
    public static class Glicko2
    {
        private const double ScaleFactor = 173.7178;
        private const double Tau = 0.5; // System constant, change as needed
        private const double Epsilon = 0.000001; // Convergence tolerance

        /// <summary>
        /// An algorithm to calculate Glicko 2 ratings adjustments.
        /// </summary>
        /// <param name="ratingA">Player A rating before recalculations</param>
        /// <param name="ratingB">Player B rating before recalculations</param>
        /// <param name="rdA">rdA should be around 350 starting, these are K values that are lowered as rating certainty is met, a min value of around 75.0 is recommended, this will result in around 16 rating lost/gained in an even matchup.</param>
        /// <param name="rdB">rdB should be around 350 starting, these are K values that are lowered as rating certainty is met, a min value of around 75.0 is recommended, this will result in around 16 rating lost/gained in an even matchup.</param>
        /// <param name="outcomeA">Player A's result: 1.0 for win 0.5 for draw 0.0 for loss</param>
        /// <returns></returns>
        public static (double newRatingA, double newRatingB) CalculateNewRatings(double ratingA, double ratingB, double rdA, double rdB, double outcomeA)
        {
            double newRatingA = Glicko2.CalculateNewRatingsA(ratingA, ratingB, rdA, rdA, outcomeA);
            double outcomeB = 0.5;
            if (outcomeA == 0.0)
                outcomeB = 1.0;
            if (outcomeA == 1.0)
                outcomeB = 0.0;
            double newRatingB = Glicko2.CalculateNewRatingsA(ratingB, ratingB, rdB, rdB, outcomeB);
            return (newRatingA, newRatingB);
        }

        public static double CalculateNewRatingsA(double ratingA, double ratingB, double rdA, double rdB, double outcome)
        {
            double muA = (ratingA - 1500) / ScaleFactor;
            double muB = (ratingB - 1500) / ScaleFactor;
            double phiA = rdA / ScaleFactor;
            double phiB = rdB / ScaleFactor;

            double gPhiB = G(phiB);
            double eA = E(muA, muB, gPhiB);

            double v = 1.0 / (gPhiB * gPhiB * eA * (1 - eA));
            double delta = v * gPhiB * (outcome - eA);

            double a = Math.Log(phiA * phiA);
            double A = a;
            double B;

            if (delta * delta > phiA * phiA + v)
            {
                B = Math.Log(delta * delta - phiA * phiA - v);
            }
            else
            {
                double k = 1;
                while (F(a - k * Tau, delta, phiA, v, a) < 0)
                {
                    k++;
                }
                B = a - k * Tau;
            }

            double fA = F(A, delta, phiA, v, a);
            double fB = F(B, delta, phiA, v, a);

            while (Math.Abs(B - A) > Epsilon)
            {
                double C = A + (A - B) * fA / (fB - fA);
                double fC = F(C, delta, phiA, v, a);

                if (fC * fB < 0)
                {
                    A = B;
                    fA = fB;
                }
                else
                {
                    fA /= 2.0;
                }

                B = C;
                fB = fC;
            }

            double newPhiA = Math.Exp(A / 2);
            double newMuA = muA + (newPhiA * newPhiA * gPhiB * (outcome - eA));

            double newRatingA = newMuA * ScaleFactor + 1500;
            double newRdA = newPhiA * ScaleFactor;

            return newRatingA; // Here we update only Player A's rating. You can similarly update Player B.
        }

        private static double G(double phi)
        {
            return 1.0 / Math.Sqrt(1.0 + 3.0 * phi * phi / (Math.PI * Math.PI));
        }

        private static double E(double muA, double muB, double gPhiB)
        {
            return 1.0 / (1.0 + Math.Exp(-gPhiB * (muA - muB)));
        }

        private static double F(double x, double delta, double phi, double v, double a)
        {
            double expX = Math.Exp(x);
            return expX * (delta * delta - phi * phi - v - expX) / (2.0 * (phi * phi + v + expX) * (phi * phi + v + expX)) - (x - a) / (Tau * Tau);
        }
    }
}
