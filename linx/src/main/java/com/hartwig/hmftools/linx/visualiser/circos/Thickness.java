package com.hartwig.hmftools.linx.visualiser.circos;

import java.util.List;

import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.linx.visualiser.data.Connector;

import org.jetbrains.annotations.NotNull;

public class Thickness
{

    private static final double MAX_TOTAL_PLOIDY = 1000;

    private final double minPixels;
    private final double maxPixels;
    private final double gradient;

    public Thickness(final double minPixels, final double maxPixels, @NotNull final List<Connector> connectors)
    {
        this.minPixels = minPixels;
        this.maxPixels = maxPixels;

        double maxPloidy = bound(connectors.stream().mapToDouble(Connector::ploidy).max().orElse(0), 6, 60);
        double unadjustedGradient = (maxPixels - minPixels) / (maxPloidy - 1);

        double totalPloidy = connectors.stream().mapToDouble(x -> thicknessPixels(x.ploidy(), unadjustedGradient)).sum();
        gradient = Doubles.greaterThan(totalPloidy, MAX_TOTAL_PLOIDY) ?
                (MAX_TOTAL_PLOIDY - connectors.size()) / (totalPloidy - connectors.size()) * unadjustedGradient
                : unadjustedGradient;

    }

    public double thicknessPixels(double ploidy)
    {
        return thicknessPixels(ploidy, gradient);
    }

    private double thicknessPixels(double ploidy, double gradient)
    {
        return Math.round(bound(gradient * ploidy + (minPixels - gradient), minPixels, maxPixels));
    }

    private static double bound(double value, double min, double max)
    {
        return Math.min(Math.max(value, min), max);
    }
}
