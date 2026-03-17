package com.hartwig.hmftools.amber.purity;

public interface GnomadFrequencySupplier
{
    double getFrequency(String chromosome, int position);
}