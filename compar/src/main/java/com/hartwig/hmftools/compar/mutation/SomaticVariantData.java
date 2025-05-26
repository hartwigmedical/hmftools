package com.hartwig.hmftools.compar.mutation;

import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_AF;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_VARIANT_CN;
import static com.hartwig.hmftools.compar.common.CommonUtils.createMismatchFromDiffs;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.variant.CodingEffect.UNDEFINED;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.PASS;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.PURPLE_BIALLELIC_FLAG;
import static com.hartwig.hmftools.common.variant.PurpleVcfTags.SUBCLONAL_LIKELIHOOD_FLAG;
import static com.hartwig.hmftools.common.variant.SageVcfTags.LOCAL_PHASE_SET;
import static com.hartwig.hmftools.common.variant.CommonVcfTags.REPORTED_FLAG;
import static com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser.VAR_IMPACT;
import static com.hartwig.hmftools.compar.common.Category.SOMATIC_VARIANT;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_QUAL;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_REPORTED;
import static com.hartwig.hmftools.compar.common.DiffFunctions.FILTER_DIFF;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkDiff;
import static com.hartwig.hmftools.compar.common.DiffFunctions.checkFilterDiffs;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_BIALLELIC;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_CANON_EFFECT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_CODING_EFFECT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_GENE;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_HGVS_CODING;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_HGVS_PROTEIN;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_HOTSPOT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_OTHER_REPORTED;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_PURITY_ADJUSTED_VAF;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_TIER;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_TUMOR_SUPPORTING_READ_COUNT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_TUMOR_TOTAL_READ_COUNT;
import static com.hartwig.hmftools.compar.mutation.VariantCommon.FLD_VARIANT_COPY_NUMBER;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.Hotspot;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;

import org.jooq.Record;

import htsjdk.variant.variantcontext.VariantContext;

public class SomaticVariantData implements ComparableItem
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;

    public final VariantType Type;
    public final String Gene;

    public final boolean Reported;
    public final Hotspot HotspotStatus;
    public final VariantTier Tier;
    public final boolean Biallelic;
    public final String CanonicalEffect;
    public final String CanonicalCodingEffect;
    public final String CanonicalHgvsCodingImpact;
    public final String CanonicalHgvsProteinImpact;
    public final String OtherReportedEffects;
    public final boolean HasLPS;
    public final int Qual;
    public final double SubclonalLikelihood;
    public final Set<String> Filters;
    public final double VariantCopyNumber;
    public final double PurityAdjustedVaf;
    public final AllelicDepth TumorDepth;
    public final boolean IsFromUnfilteredVcf;

    public String mComparisonChromosome;
    public int mComparisonPosition;

    protected static final String FLD_SUBCLONAL_LIKELIHOOD = "SubclonalLikelihood";
    protected static final String FLD_LPS = "HasLPS";

    protected static final double NO_QUAL_PRESENT = -10;

    public SomaticVariantData(
            final String chromosome, final int position, final String ref, final String alt, final VariantType type,
            final String gene, final boolean reported, final Hotspot hotspotStatus, final VariantTier tier, final boolean biallelic,
            final String canonicalEffect, final String canonicalCodingEffect, final String canonicalHgvsCodingImpact,
            final String canonicalHgvsProteinImpact, final String otherReportedEffects, final boolean hasLPS, final int qual,
            final double subclonalLikelihood, final Set<String> filters, final double variantCopyNumber, final double purityAdjustedVaf,
            final AllelicDepth tumorDepth, final boolean isFromUnfilteredVcf, final String comparisonChromosome,
            final int comparisonPosition)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        Type = type;
        Gene = gene;
        Reported = reported;
        HotspotStatus = hotspotStatus;
        Tier = tier;
        Biallelic = biallelic;
        CanonicalEffect = canonicalEffect;
        CanonicalCodingEffect = canonicalCodingEffect;
        CanonicalHgvsCodingImpact = canonicalHgvsCodingImpact;
        CanonicalHgvsProteinImpact = canonicalHgvsProteinImpact;
        OtherReportedEffects = otherReportedEffects != null ? otherReportedEffects : "";
        HasLPS = hasLPS;
        Qual = qual;
        SubclonalLikelihood = subclonalLikelihood;
        Filters = filters;
        VariantCopyNumber = variantCopyNumber;
        PurityAdjustedVaf = purityAdjustedVaf;
        TumorDepth = tumorDepth;
        IsFromUnfilteredVcf = isFromUnfilteredVcf;
        mComparisonChromosome = comparisonChromosome;
        mComparisonPosition = comparisonPosition;
    }

    public SomaticVariantData(
            final String chromosome, final int position, final String ref, final String alt, final VariantType type,
            final String gene, final boolean reported, final Hotspot hotspotStatus, final VariantTier tier, final boolean biallelic,
            final String canonicalEffect, final String canonicalCodingEffect, final String canonicalHgvsCodingImpact,
            final String canonicalHgvsProteinImpact, final String otherReportedEffects, final boolean hasLPS, final int qual,
            final double subclonalLikelihood, final Set<String> filters, final double variantCopyNumber, final double purityAdjustedVaf,
            final AllelicDepth tumorDepth, final boolean isFromUnfilteredVcf)
    {
        this(chromosome, position, ref, alt, type, gene, reported, hotspotStatus, tier, biallelic, canonicalEffect, canonicalCodingEffect,
                canonicalHgvsCodingImpact, canonicalHgvsProteinImpact, otherReportedEffects, hasLPS, qual, subclonalLikelihood, filters,
                variantCopyNumber, purityAdjustedVaf, tumorDepth, isFromUnfilteredVcf, chromosome, position);
    }

    @Override
    public Category category() { return SOMATIC_VARIANT; }

    @Override
    public String key()
    {
        if(mComparisonPosition != Position)
            return format("%s:%d %s>%s %s liftover(%s:%d)", Chromosome, Position, Ref, Alt, Type, mComparisonChromosome, mComparisonPosition);
        else
            return format("%s:%d %s>%s %s", Chromosome, Position, Ref, Alt, Type);
    }

    @Override
    public List<String> displayValues()
    {
        List<String> values = Lists.newArrayList();
        values.add(format("%s", Reported));
        values.add(format("%s", HotspotStatus));
        values.add(format("%s", Tier));
        values.add(format("%s", Biallelic));
        values.add(format("%s", Gene));
        values.add(format("%s", CanonicalEffect));
        values.add(format("%s", CanonicalCodingEffect));
        values.add(format("%s", CanonicalHgvsCodingImpact));
        values.add(format("%s", CanonicalHgvsProteinImpact));
        values.add(format("%s", OtherReportedEffects));
        values.add(format("%d", Qual));
        values.add(format("%.2f", VariantCopyNumber));
        values.add(format("%.2f", PurityAdjustedVaf));
        values.add(String.format("%d", TumorDepth.AlleleReadCount));
        values.add(String.format("%d", TumorDepth.TotalReadCount));

        values.add(format("%.2f", SubclonalLikelihood));
        values.add(format("%s", HasLPS));

        return values;
    }

    @Override
    public boolean reportable() {
        return !IsFromUnfilteredVcf && Reported;
    }

    @Override
    public boolean isPass() {
        // A reportable variant not in a gene should be impossible, but if it happens we want to see it
        return !IsFromUnfilteredVcf && (Reported || !Gene.isEmpty());
    }

    @Override
    public boolean matches(final ComparableItem other)
    {
        final SomaticVariantData otherVar = (SomaticVariantData) other;

        if(!mComparisonChromosome.equals(otherVar.Chromosome) || mComparisonPosition != otherVar.Position)
            return false;

        if(!Ref.equals(otherVar.Ref) || !Alt.equals(otherVar.Alt))
            return false;

        if(Type != otherVar.Type)
            return false;

        return true;
    }

    public String comparisonChromosome() { return mComparisonChromosome; }
    public int comparisonPosition() { return mComparisonPosition; }

    public void setComparisonCoordinates(final String chromosome, final int position)
    {
        mComparisonChromosome = chromosome;
        mComparisonPosition = position;
    }

    @Override
    public Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches)
    {
        return findMismatch(other, matchLevel, thresholds, includeMatches, false);
    }

    protected Mismatch findMismatch(final ComparableItem other, final MatchLevel matchLevel, final DiffThresholds thresholds,
            final boolean includeMatches, final boolean nonPurpleVcfs)
    {
        final SomaticVariantData otherVar = (SomaticVariantData) other;
        final List<String> diffs = findDiffs(otherVar, thresholds, nonPurpleVcfs);
        return createMismatchFromDiffs(this, other, diffs, matchLevel, includeMatches);
    }

    private List<String> findDiffs(final SomaticVariantData otherVar, final DiffThresholds thresholds, final boolean nonPurpleVcfs)
    {
        final List<String> diffs = Lists.newArrayList();

        if(Qual != NO_QUAL_PRESENT && otherVar.Qual != NO_QUAL_PRESENT)
            checkDiff(diffs, FLD_QUAL, Qual, otherVar.Qual, thresholds);

        checkDiff(diffs, FLD_REPORTED, Reported, otherVar.Reported);
        checkDiff(diffs, FLD_TIER, Tier.toString(), otherVar.Tier.toString());
        checkDiff(diffs, FLD_TUMOR_SUPPORTING_READ_COUNT, TumorDepth.AlleleReadCount, otherVar.TumorDepth.AlleleReadCount, thresholds);
        checkDiff(diffs, FLD_TUMOR_TOTAL_READ_COUNT, TumorDepth.TotalReadCount, otherVar.TumorDepth.TotalReadCount, thresholds);

        if(canComparePaveFields(otherVar))
        {
            // assumes Pave annotated - could possibly check VCF for presence of tags
            checkDiff(diffs, FLD_GENE, Gene, otherVar.Gene);
            checkDiff(diffs, FLD_CANON_EFFECT, CanonicalEffect, otherVar.CanonicalEffect);
            checkDiff(diffs, FLD_CODING_EFFECT, CanonicalCodingEffect, otherVar.CanonicalCodingEffect);
            checkDiff(diffs, FLD_HGVS_CODING, CanonicalHgvsCodingImpact, otherVar.CanonicalHgvsCodingImpact);
            checkDiff(diffs, FLD_HGVS_PROTEIN, CanonicalHgvsProteinImpact, otherVar.CanonicalHgvsProteinImpact);
        }

        if(canComparePurpleFields(otherVar, nonPurpleVcfs))
        {
            checkDiff(diffs, FLD_HOTSPOT, HotspotStatus.toString(), otherVar.HotspotStatus.toString());
            checkDiff(diffs, FLD_BIALLELIC, Biallelic, otherVar.Biallelic);
            checkDiff(diffs, FLD_OTHER_REPORTED, OtherReportedEffects, otherVar.OtherReportedEffects);
            checkDiff(diffs, FLD_SUBCLONAL_LIKELIHOOD, SubclonalLikelihood, otherVar.SubclonalLikelihood, thresholds);
            checkDiff(diffs, FLD_VARIANT_COPY_NUMBER, VariantCopyNumber, otherVar.VariantCopyNumber, thresholds);
            checkDiff(diffs, FLD_PURITY_ADJUSTED_VAF, PurityAdjustedVaf, otherVar.PurityAdjustedVaf, thresholds);
        }

        checkDiff(diffs, FLD_LPS, HasLPS, otherVar.HasLPS);

        // compare filters
        checkFilterDiffs(Filters, otherVar.Filters, diffs);

        if(Filters.isEmpty() && otherVar.Filters.isEmpty() && !diffs.contains(FILTER_DIFF))
        {
            // if ones side is filtered, suggests was filtered downstream of Sage (eg Pave or Purple) so indicate this
            if(IsFromUnfilteredVcf && !otherVar.IsFromUnfilteredVcf)
                diffs.add(format("%s(%s/%s)", FILTER_DIFF, "FILTERED", PASS));
            else if(!IsFromUnfilteredVcf && otherVar.IsFromUnfilteredVcf)
                diffs.add(format("%s(%s/%s)", FILTER_DIFF, PASS, "FILTERED"));
        }
        return diffs;
    }

    private boolean canComparePaveFields(final SomaticVariantData otherVar)
    {
        return !IsFromUnfilteredVcf && !otherVar.IsFromUnfilteredVcf;
    }

    private boolean canComparePurpleFields(final SomaticVariantData otherVar, final boolean nonPurpleVcfs)
    {
        return !IsFromUnfilteredVcf && !otherVar.IsFromUnfilteredVcf && !nonPurpleVcfs;
    }

    public static SomaticVariantData fromContext(final VariantContext context, final String sampleId, final boolean fromUnfilteredFile)
    {
        int position = context.getStart();
        String chromosome = context.getContig();
        String ref = context.getReference().getBaseString();
        String alt = !context.getAlternateAlleles().isEmpty() ? context.getAlternateAlleles().get(0).toString() : ref;

        VariantImpact variantImpact;
        if(context.hasAttribute(VAR_IMPACT))
            variantImpact = VariantImpactSerialiser.fromVariantContext(context);
        else
            variantImpact = fromSnpEffAttributes(context);

        return new SomaticVariantData(
                chromosome, position, ref, alt, VariantType.type(context),
                variantImpact.GeneName,
                context.getAttributeAsBoolean(REPORTED_FLAG, false),
                Hotspot.fromVariant(context),
                VariantTier.fromContext(context),
                context.getAttributeAsBoolean(PURPLE_BIALLELIC_FLAG, false),
                variantImpact.CanonicalEffect,
                variantImpact.CanonicalCodingEffect.toString(),
                variantImpact.CanonicalHgvsCoding,
                variantImpact.CanonicalHgvsProtein,
                variantImpact.OtherReportableEffects,
                context.hasAttribute(LOCAL_PHASE_SET),
                (int)context.getPhredScaledQual(),
                context.getAttributeAsDouble(SUBCLONAL_LIKELIHOOD_FLAG, 0),
                context.getFilters(),
                context.getAttributeAsDouble(PURPLE_VARIANT_CN, 0),
                context.getAttributeAsDouble(PURPLE_AF, 0),
                AllelicDepth.fromGenotype(context.getGenotype(sampleId)),
                fromUnfilteredFile);
    }

    public static SomaticVariantData fromRecord(final Record record)
    {
        Set<String> filters = Arrays.stream(record.getValue(SOMATICVARIANT.FILTER).split(";", -1)).collect(Collectors.toSet());
        String localPhaseSets = record.get(SOMATICVARIANT.LOCALPHASESET);
        double qual = record.getValue(Tables.SOMATICVARIANT.QUAL);
        final AllelicDepth tumorDepth =
                new AllelicDepth(record.getValue(SOMATICVARIANT.TOTALREADCOUNT), record.getValue(SOMATICVARIANT.ALLELEREADCOUNT));

        return new SomaticVariantData(
                record.getValue(Tables.SOMATICVARIANT.CHROMOSOME),
                record.getValue(Tables.SOMATICVARIANT.POSITION),
                record.getValue(Tables.SOMATICVARIANT.REF),
                record.getValue(Tables.SOMATICVARIANT.ALT),
                VariantType.valueOf(record.getValue(SOMATICVARIANT.TYPE)),
                record.getValue(Tables.SOMATICVARIANT.GENE),
                record.getValue(SOMATICVARIANT.REPORTED).intValue() == 1,
                Hotspot.valueOf(record.getValue(SOMATICVARIANT.HOTSPOT)),
                VariantTier.fromString(record.get(SOMATICVARIANT.TIER)),
                record.getValue(SOMATICVARIANT.BIALLELIC).intValue() == 1,
                record.getValue(SOMATICVARIANT.CANONICALEFFECT),
                record.getValue(SOMATICVARIANT.CANONICALCODINGEFFECT),
                record.getValue(SOMATICVARIANT.CANONICALHGVSCODINGIMPACT),
                record.getValue(SOMATICVARIANT.CANONICALHGVSPROTEINIMPACT),
                record.getValue(SOMATICVARIANT.OTHERTRANSCRIPTEFFECTS),
                localPhaseSets != null && !localPhaseSets.isEmpty(),
                (int)qual, record.getValue(SOMATICVARIANT.SUBCLONALLIKELIHOOD),
                filters,
                record.getValue(SOMATICVARIANT.VARIANTCOPYNUMBER),
                record.getValue(SOMATICVARIANT.ADJUSTEDVAF),
                tumorDepth,
                false);
    }

    private static final String SNPEFF_WORST = "SEW";
    private static final String SNPEFF_CANONICAL = "SEC";

    private static final VariantImpact INVALID_IMPACT = new VariantImpact(
            "", "", "", UNDEFINED, "", "",
            false, "", UNDEFINED, 0);

    private static VariantImpact fromSnpEffAttributes(final VariantContext context)
    {
        if(!context.hasAttribute(SNPEFF_WORST) || !context.hasAttribute(SNPEFF_CANONICAL))
            return INVALID_IMPACT;

        final List<String> worst = context.getAttributeAsStringList(SNPEFF_WORST, "");
        final List<String> canonical = context.getAttributeAsStringList(SNPEFF_CANONICAL, "");

        String canonicalGeneName = "";
        String canonicalEffect = "";
        String canonicalTranscript = "";
        CodingEffect canonicalCodingEffect = UNDEFINED;
        String canonicalHgvsCodingImpact = "";
        String canonicalHgvsProteinImpact = "";
        boolean canonicalSpliceRegion = false;
        String otherReportableEffects = "";
        CodingEffect worstCodingEffect = UNDEFINED;
        int genesAffected = 0;

        if(worst.size() == 5)
        {
            worstCodingEffect = CodingEffect.valueOf(worst.get(3));
            genesAffected = Integer.parseInt(worst.get(4));
        }

        if(canonical.size() == 6)
        {
            canonicalGeneName = canonical.get(0);
            canonicalTranscript = canonical.get(1);
            canonicalEffect = canonical.get(2);
            canonicalCodingEffect = CodingEffect.valueOf(canonical.get(3));
            canonicalHgvsCodingImpact = canonical.get(4);
            canonicalHgvsProteinImpact = canonical.get(5);

            canonicalSpliceRegion = canonicalEffect.contains("splice");
        }

        return new VariantImpact(
                canonicalGeneName, canonicalTranscript, canonicalEffect, canonicalCodingEffect, canonicalHgvsCodingImpact,
                canonicalHgvsProteinImpact, canonicalSpliceRegion, otherReportableEffects, worstCodingEffect, genesAffected);
    }

    public String toString() { return format("%s gene(%s:%s)", key(), Gene, CanonicalCodingEffect); }

    public SomaticVariantData withChromosome(String Chromosome)
    {
        return new SomaticVariantData(Chromosome, Position, Ref, Alt, Type, Gene, Reported, HotspotStatus, Tier, Biallelic, CanonicalEffect,
                CanonicalCodingEffect, CanonicalHgvsCodingImpact, CanonicalHgvsProteinImpact, OtherReportedEffects, HasLPS, Qual,
                SubclonalLikelihood, Filters, VariantCopyNumber, PurityAdjustedVaf, TumorDepth, IsFromUnfilteredVcf, mComparisonChromosome,
                mComparisonPosition);
    }

    public SomaticVariantData withPosition(int Position)
    {
        return new SomaticVariantData(Chromosome, Position, Ref, Alt, Type, Gene, Reported, HotspotStatus, Tier, Biallelic, CanonicalEffect,
                CanonicalCodingEffect, CanonicalHgvsCodingImpact, CanonicalHgvsProteinImpact, OtherReportedEffects, HasLPS, Qual,
                SubclonalLikelihood, Filters, VariantCopyNumber, PurityAdjustedVaf, TumorDepth, IsFromUnfilteredVcf, mComparisonChromosome,
                mComparisonPosition);
    }

    public SomaticVariantData withRef(String Ref)
    {
        return new SomaticVariantData(Chromosome, Position, Ref, Alt, Type, Gene, Reported, HotspotStatus, Tier, Biallelic, CanonicalEffect,
                CanonicalCodingEffect, CanonicalHgvsCodingImpact, CanonicalHgvsProteinImpact, OtherReportedEffects, HasLPS, Qual,
                SubclonalLikelihood, Filters, VariantCopyNumber, PurityAdjustedVaf, TumorDepth, IsFromUnfilteredVcf, mComparisonChromosome,
                mComparisonPosition);
    }

    public SomaticVariantData withAlt(String Alt)
    {
        return new SomaticVariantData(Chromosome, Position, Ref, Alt, Type, Gene, Reported, HotspotStatus, Tier, Biallelic, CanonicalEffect,
                CanonicalCodingEffect, CanonicalHgvsCodingImpact, CanonicalHgvsProteinImpact, OtherReportedEffects, HasLPS, Qual,
                SubclonalLikelihood, Filters, VariantCopyNumber, PurityAdjustedVaf, TumorDepth, IsFromUnfilteredVcf, mComparisonChromosome,
                mComparisonPosition);
    }

    public SomaticVariantData withType(VariantType Type)
    {
        return new SomaticVariantData(Chromosome, Position, Ref, Alt, Type, Gene, Reported, HotspotStatus, Tier, Biallelic, CanonicalEffect,
                CanonicalCodingEffect, CanonicalHgvsCodingImpact, CanonicalHgvsProteinImpact, OtherReportedEffects, HasLPS, Qual,
                SubclonalLikelihood, Filters, VariantCopyNumber, PurityAdjustedVaf, TumorDepth, IsFromUnfilteredVcf, mComparisonChromosome,
                mComparisonPosition);
    }

    public SomaticVariantData withIsFromUnfilteredVcf(boolean IsFromUnfilteredVcf)
    {
        return new SomaticVariantData(Chromosome, Position, Ref, Alt, Type, Gene, Reported, HotspotStatus, Tier, Biallelic, CanonicalEffect,
                CanonicalCodingEffect, CanonicalHgvsCodingImpact, CanonicalHgvsProteinImpact, OtherReportedEffects, HasLPS, Qual,
                SubclonalLikelihood, Filters, VariantCopyNumber, PurityAdjustedVaf, TumorDepth, IsFromUnfilteredVcf, mComparisonChromosome,
                mComparisonPosition);
    }

    public SomaticVariantData withComparisonChromosome(String ComparisonChromosome)
    {
        return new SomaticVariantData(Chromosome, Position, Ref, Alt, Type, Gene, Reported, HotspotStatus, Tier, Biallelic, CanonicalEffect,
                CanonicalCodingEffect, CanonicalHgvsCodingImpact, CanonicalHgvsProteinImpact, OtherReportedEffects, HasLPS, Qual,
                SubclonalLikelihood, Filters, VariantCopyNumber, PurityAdjustedVaf, TumorDepth, IsFromUnfilteredVcf, ComparisonChromosome,
                mComparisonPosition);
    }

    public SomaticVariantData withComparisonPosition(int ComparisonPosition)
    {
        return new SomaticVariantData(Chromosome, Position, Ref, Alt, Type, Gene, Reported, HotspotStatus, Tier, Biallelic, CanonicalEffect,
                CanonicalCodingEffect, CanonicalHgvsCodingImpact, CanonicalHgvsProteinImpact, OtherReportedEffects, HasLPS, Qual,
                SubclonalLikelihood, Filters, VariantCopyNumber, PurityAdjustedVaf, TumorDepth, IsFromUnfilteredVcf, mComparisonChromosome,
                ComparisonPosition);
    }
}
