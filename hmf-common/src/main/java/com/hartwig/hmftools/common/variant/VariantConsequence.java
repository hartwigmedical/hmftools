package com.hartwig.hmftools.common.variant;

import static com.hartwig.hmftools.common.variant.ConsequenceEffects.FIVE_PRIME_UTR_EFFECT;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.SPLICE_ACCEPTOR_EFFECT;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.SPLICE_DONOR_EFFECT;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.SPLICE_REGION_EFFECT;
import static com.hartwig.hmftools.common.variant.ConsequenceEffects.THREE_PRIME_UTR_EFFECT;

import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public enum VariantConsequence
{
    // See also http://sequenceontology.org
    // starting with the ones used by HMF
    SPLICE_ACCEPTOR_VARIANT(SPLICE_ACCEPTOR_EFFECT),
    SPLICE_DONOR_VARIANT(SPLICE_DONOR_EFFECT),
    SPLICE_REGION_VARIANT(SPLICE_REGION_EFFECT,
            "exonic_splice_region_variant", "non_coding_transcript_splice_region_variant"),
    STOP_GAINED("stop_gained"),
    STOP_LOST("stop_lost"),
    START_LOST("start_lost"),
    FRAMESHIFT_VARIANT("frameshift_variant",
            "frame_restoring_variant", "frameshift_elongation", "frameshift_truncation", "minus_1_frameshift_variant",
            "minus_2_frameshift_variant", "plus_1_frameshift_variant", "plus_2_frameshift_variant"),

    INFRAME_INSERTION("inframe_insertion", "conservative_inframe_insertion", "disruptive_inframe_insertion"),
    INFRAME_DELETION("inframe_deletion", "conservative_inframe_deletion", "disruptive_inframe_deletion"),
    MISSENSE_VARIANT("missense_variant",
            "conservative_missense_variant", "non_conservative_missense_variant",
            "rare_amino_acid_variant", "pyrrolysine_loss", "selenocysteine_loss"),
    SYNONYMOUS_VARIANT("synonymous_variant", "stop_retained_variant"),
    INTRON_VARIANT("intron_variant"),
    UTR_VARIANT("UTR_variant", THREE_PRIME_UTR_EFFECT, FIVE_PRIME_UTR_EFFECT,
            "5_prime_UTR_premature_start_codon_gain_variant", "5_prime_UTR_truncation", "3_prime_UTR_truncation"),
    UPSTREAM_GENE_VARIANT("upstream_gene_variant"),

    INTRAGENIC_VARIANT("intragenic_variant"),
    TRANSCRIPT("transcript"),
    NON_CODING_TRANSCRIPT_VARIANT("non_coding_transcript_variant", "non_coding_transcript_exon_variant"),
    SEQUENCE_FEATURE("sequence_feature"),
    REGULATORY_REGION_VARIANT("regulatory_region_variant", "TF_binding_site_variant"),
    INITIATOR_CODON_VARIANT("initiator_codon_variant"),
    EXON_LOSS_VARIANT("exon_loss_variant", "exon_loss"),
    NON_CANONICAL_START_CODON("non_canonical_start_codon"),
    TRANSCRIPT_ABLATION("transcript_ablation"),
    TFBS_ABLATION("TFBS_ablation"),
    TRANSCRIPT_AMPLIFICATION("transcript_amplification"),
    STRUCTURAL_INTERACTION_VARIANT("structural_interaction_variant"),
    FUSION("bidirectional_gene_fusion", "gene_fusion"),
    PROTEIN_PROTEIN_CONTACT("protein_protein_contact"),
    OTHER(Strings.EMPTY);


    private final String mParentSequenceOntologyTerm;
    private final List<String> mSequenceOntologySubTerms;

    VariantConsequence(final String parentSequenceOntologyTerm, final String... sequenceOntologySubTerms)
    {
        mParentSequenceOntologyTerm = parentSequenceOntologyTerm;
        mSequenceOntologySubTerms = Lists.newArrayList(sequenceOntologySubTerms);
    }

    public static List<VariantConsequence> convertFromEffects(@NotNull final List<String> effects)
    {
        final List<VariantConsequence> consequences = Lists.newArrayList();

        effects.forEach(x -> consequences.add(fromEffect(x)));
        return consequences;
    }

    public static VariantConsequence fromEffect(@NotNull final String effect)
    {
        for(final VariantConsequence consequence : VariantConsequence.values())
        {
            if(consequence.isParentTypeOf(effect))
                return consequence;
        }

        return VariantConsequence.OTHER;
    }

    public static String consequencesToString(final List<VariantConsequence> consequences, final String delim)
    {
        StringJoiner sj = new StringJoiner(delim);
        consequences.forEach(x -> sj.add(x.toString()));
        return sj.toString();
    }

    public static String consequenceString(final List<VariantConsequence> consequences)
    {
        final StringJoiner consequenceString = new StringJoiner("; ");
        for(final VariantConsequence consequence : consequences)
        {
            if(!consequence.readableSequenceOntologyTerm().isEmpty())
            {
                consequenceString.add(consequence.readableSequenceOntologyTerm());
            }
        }
        return consequenceString.toString();
    }

    public boolean isParentTypeOf(@NotNull final String annotation)
    {
        return annotation.equals(mParentSequenceOntologyTerm) || mSequenceOntologySubTerms.contains(annotation);
    }

    public String description()
    {
        if(!mSequenceOntologySubTerms.isEmpty())
            return mSequenceOntologySubTerms.get(0);

        return mParentSequenceOntologyTerm;
    }

    public int rank()
    {
        switch(this)
        {
            case STOP_LOST:
            case STOP_GAINED:
            case START_LOST:
                return 60;

            case SPLICE_DONOR_VARIANT:
            case SPLICE_ACCEPTOR_VARIANT:
                return 50;

            case SPLICE_REGION_VARIANT: return 49;

            case FRAMESHIFT_VARIANT: return 40;
            case MISSENSE_VARIANT: return 40;

            case INFRAME_INSERTION:
            case INFRAME_DELETION:
                return 30;

            case SYNONYMOUS_VARIANT: return 25;

            case UTR_VARIANT: return 20;
            case INTRON_VARIANT: return 20;

            case UPSTREAM_GENE_VARIANT: return 10;

            default:
                return 0;
        }
    }

    @NotNull
    public String readableSequenceOntologyTerm()
    {
        return mParentSequenceOntologyTerm.replace("_", " ");
    }
}
