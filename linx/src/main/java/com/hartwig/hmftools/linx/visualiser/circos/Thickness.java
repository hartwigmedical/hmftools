package com.hartwig.hmftools.linx.visualiser.circos;

import java.util.List;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.linx.visualiser.data.Connector;

import org.jetbrains.annotations.NotNull;

public class Thickness
{
    private static final double MAX_TOTAL_PLOIDY = 1000;

    private final double mMinPixels;
    private final double mMaxPixels;
    private final double mGradient;

    public Thickness(final double minPixels, final double maxPixels, @NotNull final List<Connector> connectors)
    {
        mMinPixels = minPixels;
        mMaxPixels = maxPixels;

        double maxPloidy = bound(connectors.stream().mapToDouble(Connector::ploidy).max().orElse(0), 6, 60);
        double unadjustedGradient = (maxPixels - minPixels) / (maxPloidy - 1);

        double totalPloidy = connectors.stream().mapToDouble(x -> thicknessPixels(x.ploidy(), unadjustedGradient)).sum();
        mGradient = Doubles.greaterThan(totalPloidy, MAX_TOTAL_PLOIDY) ?
                (MAX_TOTAL_PLOIDY - connectors.size()) / (totalPloidy - connectors.size()) * unadjustedGradient
                : unadjustedGradient;

    }

    public double thicknessPixels(double ploidy)
    {
        return thicknessPixels(ploidy, mGradient);
    }

    private double thicknessPixels(double ploidy, double gradient)
    {
        return Math.round(bound(gradient * ploidy + (mMinPixels - gradient), mMinPixels, mMaxPixels));
    }

    private static double bound(double value, double min, double max)
    {
        return Math.min(Math.max(value, min), max);
    }
}
