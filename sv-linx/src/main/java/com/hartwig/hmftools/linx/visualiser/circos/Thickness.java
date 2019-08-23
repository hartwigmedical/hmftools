package com.hartwig.hmftools.linx.visualiser.circos;

import java.util.List;

import com.hartwig.hmftools.common.numeric.Doubles;
import com.hartwig.hmftools.linx.visualiser.data.Connector;

import org.jetbrains.annotations.NotNull;

public class Thickness
{

    private static final double MAX_TOTAL_PLOIDY = 1000;

    private static final double MIN_PIXELS = 1;
    private static final double MAX_PIXELS = 12;

    private final double gradient;

    public Thickness(@NotNull final List<Connector> connectors)
    {
        double maxPloidy = bound(connectors.stream().mapToDouble(Connector::ploidy).max().orElse(0), 6, 60);

        double maxPixelGradient = (MAX_PIXELS - MIN_PIXELS) / (maxPloidy - 1);
        double totalPloidy = connectors.stream().mapToDouble(x -> thicknessPixels(x.ploidy(), maxPixelGradient)).sum();
        gradient = Doubles.greaterThan(totalPloidy, MAX_TOTAL_PLOIDY) ?
                (MAX_TOTAL_PLOIDY - connectors.size()) / (totalPloidy - connectors.size()) * maxPixelGradient
                : maxPixelGradient;

    }

    public double thicknessPixels(double ploidy)
    {
        return thicknessPixels(ploidy, gradient);
    }

    private static double thicknessPixels(double ploidy, double gradient)
    {
        return Math.round(bound(gradient * ploidy + (1 - gradient), MIN_PIXELS, MAX_PIXELS));
    }

    private static double bound(double value, double min, double max)
    {
        return Math.min(Math.max(value, min), max);
    }
}
