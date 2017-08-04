package com.hartwig.hmftools.common.variant.strelka;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.exception.HartwigException;
import com.hartwig.hmftools.common.variant.Variant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.apache.logging.log4j.util.Strings;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class StrelkaSomaticVariant implements Variant {
    private static final String SNP_QUAL_FIELD = "QSS_NT";
    private static final String SNP_TIER_INDEX_FIELD = "TQSS_NT";
    private static final String INDEL_QUAL_FIELD = "QSI_NT";
    private static final String INDEL_TIER_INDEX_FIELD = "TQSI_NT";
    private static final String TIR_FIELD = "TIR";
    private static final String TAR_FIELD = "TAR";
    private static final String ALLELE_SEPARATOR = ",";

    @NotNull
    public abstract String originalVCFLine();

    @Override
    @NotNull
    public abstract VariantType type();

    @Override
    @NotNull
    public abstract String chromosome();

    @Override
    public abstract long position();

    @Override
    @NotNull
    public abstract String ref();

    @Override
    @NotNull
    public abstract String alt();

    @Override
    @NotNull
    public abstract String filter();

    @NotNull
    public abstract Map<String, String> infoData();

    @NotNull
    public abstract ArrayListMultimap<String, String> tumorData();

    private int tierIndex() throws HartwigException {
        if (type() == VariantType.SNP) {
            return Integer.parseInt(infoData().get(SNP_TIER_INDEX_FIELD)) - 1;
        } else if (type() == VariantType.INDEL) {
            return Integer.parseInt(infoData().get(INDEL_TIER_INDEX_FIELD)) - 1;
        } else {
            throw new HartwigException("record is not indel or snp: " + originalVCFLine());
        }
    }

    public int qualityScore() throws HartwigException {
        if (type() == VariantType.SNP) {
            return Integer.parseInt(infoData().get(SNP_QUAL_FIELD));
        } else if (type() == VariantType.INDEL) {
            return Integer.parseInt(infoData().get(INDEL_QUAL_FIELD));
        } else {
            throw new HartwigException("record is not indel or snp: " + originalVCFLine());
        }
    }

    public double allelicFrequency() throws HartwigException {
        if (type() == VariantType.SNP) {
            return readAf(alt() + "U", ref() + "U");
        } else if (type() == VariantType.INDEL) {
            return readAf(TIR_FIELD, TAR_FIELD);
        } else {
            throw new HartwigException("record is not indel or snp: " + originalVCFLine());
        }
    }

    private double readAf(@NotNull final String first, @NotNull final String second) throws HartwigException {
        final double firstField = Double.parseDouble(tumorData().get(first).get(tierIndex()));
        final double secondField = Double.parseDouble(tumorData().get(second).get(tierIndex()));
        final double sum = firstField + secondField;
        if (sum == 0) {
            return 0;
        }
        return firstField / sum;
    }

    public String readRefAD() {
        if (type() == VariantType.SNP) {
            return readAD(ref());
        } else if (type() == VariantType.INDEL) {
            return tumorData().get(TAR_FIELD).get(0);
        } else {
            return "0";
        }
    }

    public String readAltAD() {
        if (type() == VariantType.SNP) {
            return readAD(alt());
        } else if (type() == VariantType.INDEL) {
            return tumorData().get(TIR_FIELD).get(0);
        } else {
            return "0";
        }
    }

    @NotNull
    private String readAD(@NotNull final String allelesField) {
        final List<String> allAlleles = Lists.newArrayList(allelesField.split(ALLELE_SEPARATOR));
        final List<String> alleleAds = allAlleles.stream().map(allele -> {
            final String alleleKey = allele + "U";
            if (tumorData().get(alleleKey).size() > 0) {
                final String alleleAd = tumorData().get(alleleKey).get(0);
                if (!alleleAd.equals("") && !alleleAd.equals(".")) {
                    return alleleAd;
                }
            }
            return "0";
        }).collect(Collectors.toList());
        return Strings.join(alleleAds, ',');
    }
}
