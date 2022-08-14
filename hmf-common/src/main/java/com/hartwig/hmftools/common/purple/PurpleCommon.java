package com.hartwig.hmftools.common.purple;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

import org.jetbrains.annotations.NotNull;

public final class PurpleCommon
{
    public static final String PURPLE_SOMATIC_VCF_SUFFIX = ".purple.somatic.vcf.gz";
    public static final String PURPLE_GERMLINE_VCF_SUFFIX = ".purple.germline.vcf.gz";
    public static final String PURPLE_SV_VCF_SUFFIX = ".purple.sv.vcf.gz";

    public static final String DELIMITER = "\t";

    public static DecimalFormat decimalFormat(@NotNull String format)
    {
        return new DecimalFormat(format, DecimalFormatSymbols.getInstance(Locale.ENGLISH));
    }
}
