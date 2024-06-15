package com.hartwig.hmftools.cup.common;

import static com.hartwig.hmftools.common.variant.VariantType.INDEL;
import static com.hartwig.hmftools.common.variant.VariantType.SNP;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.cup.feature.KnownMutation;

public class CupConstants
{
    public static final String APP_NAME = "Cuppa";

    public static final int GEN_POS_BUCKET_SIZE = 500000;
    public static final int GEN_POS_MAX_SAMPLE_COUNT = 20000;

    // common
    public static final double MIN_CLASSIFIER_SCORE = 0.01;
    public static final double FEATURE_DAMPEN_FACTOR_DEFAULT = 0.75;

    public static final String CANCER_TYPE_UNKNOWN = "Unknown";
    public static final String CANCER_TYPE_OTHER = "Other";

    public static final List<String> AID_APOBEC_TRINUCLEOTIDE_CONTEXTS = Lists.newArrayList(
            "C>T_TCA", "C>T_TCC", "C>T_TCG", "C>T_TCT", "C>G_TCA", "C>G_TCC", "C>G_TCG", "C>G_TCT");

    // cancer types with gender-exclusions
    public static final String CANCER_TYPE_PROSTATE = "Prostate";
    public static final String CANCER_TYPE_OVARY = "Ovary";
    public static final String CANCER_TYPE_UTERUS = "Uterus";
    public static final String CANCER_TYPE_TESTIS = "Testis";
    public static final String CANCER_TYPE_BREAST = "Breast";
    public static final String CANCER_TYPE_BREAST_TRIPLE_NEGATIVE = "Breast triple negative";

    public static final List<KnownMutation> KNOWN_MUTATIONS = Lists.newArrayList();

    public static void loadKnownMutations(final RefGenomeVersion refGenomeVersion)
    {
        if(refGenomeVersion.is37())
        {
            // p.Thr790Met
            KNOWN_MUTATIONS.add(new KnownMutation("EGFR", SNP, "C", "T", 55249071, 55249071));

            // p.Leu858Ar
            KNOWN_MUTATIONS.add(new KnownMutation("EGFR", SNP, "T", "G", 55259515, 55259515));

            // inframe DEL in exon 19 (canonical transcript)
            KNOWN_MUTATIONS.add(new KnownMutation("EGFR", INDEL, "", "", 55242415, 55242513));

            // exon 20
            KNOWN_MUTATIONS.add(new KnownMutation("EGFR", INDEL, "", "", 55248986, 55249171));
        }
        else
        {
            KNOWN_MUTATIONS.add(new KnownMutation("EGFR", SNP, "C", "T", 55181378, 55181378));
            KNOWN_MUTATIONS.add(new KnownMutation("EGFR", SNP, "T", "G", 55191822, 55191822));
            KNOWN_MUTATIONS.add(new KnownMutation("EGFR", INDEL, "", "", 55174722, 55174820));
            KNOWN_MUTATIONS.add(new KnownMutation("EGFR", INDEL, "", "", 55181293, 55181478));
        }
    }

    public static boolean isCandidateCancerType(final Gender gender, final String cancerType)
    {
        if(gender == null)
            return true;

        if(cancerType.contains(CANCER_TYPE_UTERUS) || cancerType.contains(CANCER_TYPE_OVARY))
        {
            if(gender != Gender.FEMALE)
                return false;
        }

        if(cancerType.contains(CANCER_TYPE_PROSTATE) || cancerType.contains(CANCER_TYPE_TESTIS))
        {
            if(gender == Gender.FEMALE)
                return false;
        }

        return true;
    }
}
