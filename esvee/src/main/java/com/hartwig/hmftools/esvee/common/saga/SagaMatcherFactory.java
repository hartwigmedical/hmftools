package com.hartwig.hmftools.esvee.common.saga;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.util.Objects.requireNonNull;
import static java.util.function.UnaryOperator.identity;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.bwa.BwaMemAligner;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

import htsjdk.samtools.SAMSequenceRecord;

// For efficient creation of SagaMatcher instances.
public class SagaMatcherFactory
{
    private final SagaResource mSagaResource;

    private final SagaLocationMatcher.Config mLocationMatcherConfig;
    private final Map<String, List<SagaIndexedBreakend>> mSearchableBreakends;

    private final SagaSequenceMatcher.Config mSequenceMatcherConfig;
    private final Map<Integer, SagaAssembly> mAssembliesByContigId;

    public SagaMatcherFactory(final SagaResource sagaResource, final SagaLocationMatcher.Config locationMatcherConfig,
            final SagaSequenceMatcher.Config sequenceMatcherConfig)
    {
        mSagaResource = sagaResource;

        mSearchableBreakends = sagaResource.assemblies().stream()
                .flatMap(assembly ->
                        assembly.variant().breakends().map(breakend ->
                                new SagaIndexedBreakend(breakend, assembly.variant())))
                .collect(Collectors.groupingBy(SagaIndexedBreakend::chromosome));
        // Sort by position so they can be binary-searched.
        mSearchableBreakends.forEach((chr, breakends) -> breakends.sort(null));

        mLocationMatcherConfig = locationMatcherConfig;

        mSequenceMatcherConfig = sequenceMatcherConfig;
        // This is precomputed to avoid creating it once per instance.
        Map<Integer, String> contigIdToName = sagaResource.samDict().getSequences()
                .stream()
                .collect(Collectors.toMap(SAMSequenceRecord::getSequenceIndex, SAMSequenceRecord::getSequenceName));
        Map<String, SagaAssembly> contigToAssembly = sagaResource.assemblies().stream()
                .collect(Collectors.toMap(SagaAssembly::fastaLabel, identity()));
        mAssembliesByContigId = contigIdToName.entrySet().stream()
                .collect(Collectors.toMap(Map.Entry::getKey, entry -> requireNonNull(contigToAssembly.get(entry.getValue()))));
    }

    public SagaMatcherFactory(final String sagaFastaFile)
    {
        this(new SagaResource(sagaFastaFile));
    }

    public SagaMatcherFactory(final SagaResource sagaResource)
    {
        this(sagaResource, SagaLocationMatcher.Config.DEFAULT, SagaSequenceMatcher.Config.DEFAULT);
    }

    // Creates a SagaLocationMatcher instance for matching variants within the specified region only.
    public SagaLocationMatcher createLocationMatcher(final ChrBaseRegion subregion)
    {
        return new SagaLocationMatcher(mLocationMatcherConfig, sliceBreakends(subregion));
    }

    private Map<String, List<SagaIndexedBreakend>> sliceBreakends(ChrBaseRegion subregion)
    {
        // Need to find all variants within locationDistanceMax bases, so need to expand the region by that much.
        subregion = new ChrBaseRegion(subregion.chromosome(),
                subregion.start() - mLocationMatcherConfig.locationDistanceMax(),
                subregion.end() + mLocationMatcherConfig.locationDistanceMax());

        List<SagaIndexedBreakend> chrBreakends = mSearchableBreakends.get(subregion.chromosome());
        if(chrBreakends == null || chrBreakends.isEmpty())
        {
            return Map.of();
        }

        int startIndex = max(SagaIndexedBreakend.binarySearch(chrBreakends, subregion.start()).left - 1, 0);
        int endIndex = min(SagaIndexedBreakend.binarySearch(chrBreakends, subregion.end()).right + 1, chrBreakends.size());

        chrBreakends = chrBreakends.subList(startIndex, endIndex);
        return Map.of(subregion.chromosome(), chrBreakends);
    }

    public SagaSequenceMatcher createSequenceMatcher()
    {
        // Must create one aligner per matcher, particularly because I don't think the aligner instance is thread safe.
        BwaMemAligner aligner =
                new BwaMemAligner(mSagaResource.bwaIndexImagePath(), SagaSequenceMatcher.createAlignerParams(mSequenceMatcherConfig));
        return new SagaSequenceMatcher(mSequenceMatcherConfig, aligner, mAssembliesByContigId);
    }
}
