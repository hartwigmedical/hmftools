package com.hartwig.hmftools.common.purple;

import java.util.Arrays;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

public enum PurpleQCStatus
{
    PASS,

    WARN_DELETED_GENES,
    WARN_HIGH_COPY_NUMBER_NOISE,
    WARN_GENDER_MISMATCH,
    WARN_LOW_PURITY,

    FAIL_CONTAMINATION,
    FAIL_NO_TUMOR;

    // default QC status determination
    public static final int MAX_UNSUPPORTED_SEGMENTS = 220;
    public static final int MAX_DELETED_GENES = 280;
    public static final double MAX_CONTAMINATION = 0.1;
    public static final double MIN_PURITY = 0.2;

    public static final String STATUS_DELIM = ",";

    public static boolean genderPass(final Gender amberGender, final Gender cobaltGender, final Set<GermlineAberration> aberrations)
    {
        return cobaltGender == amberGender || aberrations.contains(GermlineAberration.KLINEFELTER);
    }

    public static Set<PurpleQCStatus> calcStatus(
            boolean genderPass, int unsupportedCopyNumberSegments, int deletedGenes, double purity, final FittedPurityMethod method,
            double contamination, final int maxDeletedGenes)
    {
        Set<PurpleQCStatus> result = Sets.newHashSet();

        if(!genderPass)
        {
            result.add(PurpleQCStatus.WARN_GENDER_MISMATCH);
        }

        if(unsupportedCopyNumberSegments > MAX_UNSUPPORTED_SEGMENTS)
        {
            result.add(PurpleQCStatus.WARN_HIGH_COPY_NUMBER_NOISE);
        }

        if(deletedGenes > maxDeletedGenes)
        {
            result.add(PurpleQCStatus.WARN_DELETED_GENES);
        }

        if(Doubles.lessThan(purity, MIN_PURITY) && !method.equals(FittedPurityMethod.NO_TUMOR))
        {
            result.add(PurpleQCStatus.WARN_LOW_PURITY);
        }

        if(Doubles.greaterThan(contamination, MAX_CONTAMINATION))
        {
            result.add(PurpleQCStatus.FAIL_CONTAMINATION);
        }

        if(method.equals(FittedPurityMethod.NO_TUMOR))
        {
            result.add(PurpleQCStatus.FAIL_NO_TUMOR);
        }

        if(result.isEmpty())
        {
            result.add(PurpleQCStatus.PASS);
        }

        return result;
    }

    @NotNull
    public static String toString(final Set<PurpleQCStatus> status)
    {
        return status.stream().map(Enum::toString).collect(Collectors.joining(STATUS_DELIM));
    }

    @NotNull
    public static Set<PurpleQCStatus> fromString(final String line)
    {
        return Arrays.stream(line.split(STATUS_DELIM)).map(PurpleQCStatus::statusFromString).collect(Collectors.toSet());
    }

    private static PurpleQCStatus statusFromString(String value)
    {
        switch(value)
        {
            case "FAIL_SEGMENT":
                return WARN_HIGH_COPY_NUMBER_NOISE;
            case "FAIL_GENDER":
                return WARN_GENDER_MISMATCH;
            case "FAIL_DELETED_GENES":
                return WARN_DELETED_GENES;
        }

        return PurpleQCStatus.valueOf(value);
    }

}