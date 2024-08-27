package com.hartwig.hmftools.purple.sv;

import static com.hartwig.hmftools.common.sv.SvFactoryInterface.buildSvFactory;
import static com.hartwig.hmftools.purple.PurpleConstants.RECOVERY_MIN_LENGTH;
import static com.hartwig.hmftools.purple.PurpleConstants.RECOVERY_MIN_MATE_UNCERTAINTY;
import static com.hartwig.hmftools.purple.PurpleConstants.RECOVERY_MIN_PLOIDY;
import static com.hartwig.hmftools.purple.PurpleConstants.RECOVERY_MIN_PLOIDY_PERC;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.gripss.GripssFilters;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.sv.SvFactoryInterface;
import com.hartwig.hmftools.common.variant.GenotypeIds;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.fitting.PurityAdjuster;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegPloidy;
import com.hartwig.hmftools.purple.copynumber.sv.StructuralVariantLegPloidyFactory;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.jetbrains.annotations.Nullable;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;

public class RecoveredVariantFactory implements AutoCloseable
{
    public static final Set<String> DO_NOT_RESCUE =
            Sets.newHashSet("af", "qual", GripssFilters.DEDUP, GripssFilters.MIN_TUMOR_AF);

    private final AbstractFeatureReader<VariantContext,LineIterator> mReader;
    private final SvFactoryInterface mSvFactory;
    private final StructuralVariantLegPloidyFactory<PurpleCopyNumber> mPloidyFactory;
    private final int mMinMateQual;
    private final int mMinSglQual;

    public RecoveredVariantFactory(
            final PurityAdjuster purityAdjuster, final String recoveryVCF, final PurpleConfig config)
    {
        mReader = getFeatureReader(recoveryVCF, new VCFCodec(), true);
        VCFHeader vcfHeader = (VCFHeader)mReader.getHeader();
        GenotypeIds genotypeIds = GenotypeIds.fromVcfHeader(vcfHeader, config.ReferenceId, config.TumorId);

        mSvFactory = buildSvFactory(config.UseGridssSVs, new CompoundFilter(false));
        mSvFactory.setGenotypeOrdinals(genotypeIds.ReferenceOrdinal, genotypeIds.TumorOrdinal);
        mPloidyFactory = new StructuralVariantLegPloidyFactory<>(purityAdjuster, PurpleCopyNumber::averageTumorCopyNumber);
        mMinMateQual = config.Fitting.RecoveryMinMateQualScore;
        mMinSglQual = config.Fitting.RecoveryMinSglQualScore;
    }

    protected List<RecoveredVariant> recoverVariantAtIndex(
            int expectedOrientation, double unexplainedCopyNumberChange, int index, final List<PurpleCopyNumber> copyNumbers) throws IOException
    {
        if(index <= 0)
            return Lists.newArrayList();

        return recoverAllVariantAtIndex(expectedOrientation, unexplainedCopyNumberChange, index, copyNumbers);
    }

    protected static Optional<RecoveredVariant> findTopCandidate(final List<RecoveredVariant> variants)
    {
        if(variants.isEmpty())
            return Optional.empty();

        // priorise by non-SGLs over SGLs, then highest qual-score
        RecoveredVariant topVariant = null;

        for(RecoveredVariant variant : variants)
        {
            if(topVariant == null)
            {
                topVariant = variant;
                continue;
            }


            else if(topVariant.mate() == null && variant.mate() != null)
                topVariant = variant;
            else if(variant.context().getPhredScaledQual() > topVariant.context().getPhredScaledQual())
                topVariant = variant;
        }

        return Optional.of(topVariant);
    }

    private List<RecoveredVariant> recoverAllVariantAtIndex(
            int expectedOrientation, double unexplainedCopyNumberChange, int index, final List<PurpleCopyNumber> copyNumbers) throws IOException
    {
        assert (index > 1);

        final List<RecoveredVariant> result = Lists.newArrayList();

        final PurpleCopyNumber prev = copyNumbers.get(index - 1);
        final PurpleCopyNumber current = copyNumbers.get(index);
        final long minPosition = Math.max(1, current.minStart() - 1000);
        final long maxPosition = current.maxStart() + 1000;

        final List<VariantContext> recovered = findVariants(current.chromosome(), minPosition, maxPosition);
        for(VariantContext potentialVariant : recovered)
        {
            final String alt = potentialVariant.getAlternateAllele(0).getDisplayString();

            final String mateId = StructuralVariantFactory.mateId(potentialVariant);
            final String mateLocation = mateLocation(alt);
            final String mateChromosome = mateChromosome(mateLocation);
            final Long matePosition = matePosition(mateLocation);
            final int uncertainty = uncertainty(potentialVariant);

            final boolean mateExists = mateChromosome != null && matePosition != null && mateId != null;
            if(mateExists && !HumanChromosome.contains(mateChromosome))
                continue;

            final VariantContext mate = mateExists
                    ? findMate(mateId, mateChromosome, Math.max(1, matePosition - uncertainty), matePosition + uncertainty)
                    : null;

            if(mate != null && !isAppropriatelyFiltered(mate))
                continue;

            final StructuralVariant sv = mate != null ?
                    mSvFactory.createSV(potentialVariant, mate) : mSvFactory.createSingleBreakend(potentialVariant);

            final Optional<StructuralVariantLegPloidy> structuralVariantLegPloidy =
                    mPloidyFactory.singleLegPloidy(sv.start(), prev.averageTumorCopyNumber(), current.averageTumorCopyNumber());

            if(structuralVariantLegPloidy.isPresent())
            {
                double ploidy = structuralVariantLegPloidy.get().averageImpliedPloidy();

                if(sv.start().orientation() == expectedOrientation
                && sufficientPloidy(ploidy, unexplainedCopyNumberChange)
                && hasPotential(sv, copyNumbers))
                {
                    result.add(ImmutableRecoveredVariant.builder()
                            .context(potentialVariant)
                            .mate(mate)
                            .variant(sv)
                            .copyNumber(current)
                            .prevCopyNumber(prev)
                            .build());
                }
            }
        }

        return result;
    }

    private boolean sufficientPloidy(double ploidy, double unexplainedCopyNumberChange)
    {
        return Doubles.greaterOrEqual(ploidy, unexplainedCopyNumberChange * RECOVERY_MIN_PLOIDY_PERC)
                && Doubles.greaterOrEqual(ploidy, RECOVERY_MIN_PLOIDY);
    }

    private int uncertainty(final VariantContext context)
    {
        final int homlen = 2 * context.getAttributeAsInt("HOMLEN", 0);
        final int cipos = cipos(context);
        return Math.max(homlen, cipos);
    }

    private int cipos(final VariantContext context)
    {
        int max = RECOVERY_MIN_MATE_UNCERTAINTY;
        if(context.hasAttribute("IMPRECISE"))
        {

            final String cipos = context.getAttributeAsString("CIPOS", "-0,0");
            if(cipos.contains(","))
            {
                for(String s : cipos.split(","))
                {
                    try
                    {
                        max = Math.max(max, 2 * Math.abs(Integer.parseInt(s)));
                    } catch (Exception ignored)
                    {

                    }
                }
            }
        }
        return max;
    }

    private boolean hasPotential(final StructuralVariant variant, final List<PurpleCopyNumber> copyNumbers)
    {
        StructuralVariantLeg start = variant.start();
        StructuralVariantLeg end = variant.end();

        // This should never actually occur because we are searching within this area
        if(!isInRangeOfCopyNumberSegment(start, copyNumbers))
            return false;

        if(end == null)
            return variant.qualityScore() >= mMinSglQual;

        if(variant.qualityScore() < mMinMateQual)
            return false;

        long endPosition = end.position();
        StructuralVariantType type = variant.type();
        if(type == StructuralVariantType.DEL || type == StructuralVariantType.DUP || type == StructuralVariantType.INS)
        {
            assert (variant.end() != null);

            long length = Math.abs(endPosition - variant.start().position());
            return length >= RECOVERY_MIN_LENGTH;
        }

        return true;
    }

    private boolean isInRangeOfCopyNumberSegment(final StructuralVariantLeg leg, final List<PurpleCopyNumber> copyNumbers)
    {
        final Predicate<PurpleCopyNumber> chrRange = copyNumber -> copyNumber.chromosome().equals(leg.chromosome());
        final Predicate<PurpleCopyNumber> posRange =
                copyNumber -> leg.cnaPosition() >= copyNumber.minStart() - 1000 && leg.cnaPosition() <= copyNumber.maxStart() + 1000;
        return copyNumbers.stream().anyMatch(chrRange.and(posRange));
    }

    private List<VariantContext> findVariants(final String chromosome, final long lowerBound, final long upperBound)
            throws IOException
    {
        return mReader.query(chromosome, (int) lowerBound, (int) upperBound)
                .stream()
                .filter(RecoveredVariantFactory::isAppropriatelyFiltered)
                .collect(Collectors.toList());
    }

    @VisibleForTesting
    static boolean isAppropriatelyFiltered(VariantContext variantContext)
    {
        final Set<String> filters = variantContext.getFilters();
        return !filters.isEmpty() && filters.stream().noneMatch(DO_NOT_RESCUE::contains);
    }

    private VariantContext findMate(final String id, final String chromosome, final long min, final long max)
            throws IOException
    {
        return mReader.query(chromosome, (int) min, (int) max)
                .stream()
                .filter(x -> x.getID().equals(id))
                .findFirst()
                .orElseThrow(() -> new IOException("Unable to find mateId " + id + " between " + min + " and " + max));
    }

    static String mateLocation(final String alt)
    {
        final String bracket;
        if(alt.contains("["))
        {
            bracket = "\\[";
        }
        else if(alt.contains("]"))
        {
            bracket = "]";
        }
        else
        {
            return null;
        }

        String[] results = alt.split(bracket);
        for(String result : results)
        {
            if(result.contains(":"))
            {
                return result;
            }
        }
        return null;
    }

    static String mateChromosome(@Nullable String mate)
    {
        return mate == null || !mate.contains(":") ? null : mate.split(":")[0];
    }

    private static Long matePosition(@Nullable String mate)
    {
        return mate == null || !mate.contains(":") ? null : Long.valueOf(mate.split(":")[1]);
    }

    @Override
    public void close() throws IOException
    {
        mReader.close();
    }
}
