package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.types.BachelorConfig.BACH_LOGGER;

import java.io.File;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotEvidence;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspotEvidenceFactory;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

class BamCountReader
{
    private IndexedFastaSequenceFile mIndexedFastaSequenceFile;
    private File mRefGenomeFile;
    private SamReader mTumorReader;

    private static final int DEFAULT_MIN_BASE_QUALITY = 13;
    private static final int DEFAULT_MIN_MAPPING_QUALITY = 1;

    BamCountReader()
    {
        mIndexedFastaSequenceFile = null;
        mTumorReader = null;
        mRefGenomeFile = null;
    }

    void initialise(final String refGenomeFile, IndexedFastaSequenceFile ifSeqFile)
    {
        mIndexedFastaSequenceFile = ifSeqFile;

        mRefGenomeFile = new File(refGenomeFile);
    }

    void readBamCounts(final String bamFile, List<BachelorGermlineVariant> bachRecords)
    {
        BACH_LOGGER.debug("Reading BAM file: {}", bamFile);

        mTumorReader = SamReaderFactory.makeDefault().referenceSequence(mRefGenomeFile).open(new File(bamFile));

        final Set<VariantHotspot> allHotspots = Sets.newHashSet();

        for(BachelorGermlineVariant variant : bachRecords)
        {
            VariantHotspot variantHotspot = ImmutableVariantHotspotImpl.builder()
                    .chromosome(variant.Chromosome)
                    .position(variant.Position)
                    .ref(variant.Ref)
                    .alt(variant.Alts)
                    .build();

            allHotspots.add(variantHotspot);
        }

        final VariantHotspotEvidenceFactory hotspotEvidenceFactory = new VariantHotspotEvidenceFactory(DEFAULT_MIN_MAPPING_QUALITY, DEFAULT_MIN_BASE_QUALITY, allHotspots);
        final List<VariantHotspotEvidence> tumorEvidence = hotspotEvidenceFactory.evidence(mIndexedFastaSequenceFile, mTumorReader);

        /*
        if(tumorEvidence.size() != bachRecords.size())
        {
            BACH_LOGGER.error("Incomplete BAM evidence read: evidenceCount({}) vs bachRecords({})", tumorEvidence.size(), bachRecords.size());
            return;
        }
        */

        for(BachelorGermlineVariant variant : bachRecords)
        {
            for(VariantHotspotEvidence evidence : tumorEvidence)
            {
                if(evidence.chromosome().equals(variant.Chromosome) && evidence.position() == variant.Position
                && evidence.ref().equals(variant.Ref) && evidence.alt().equals(variant.Alts))
                {
                    variant.setTumorData(evidence.altSupport(), evidence.readDepth());

                    BACH_LOGGER.debug("Chr({}) position({}) matched, counts(ref={} alt={} depth={})",
                            variant.Chromosome, variant.Position,
                            variant.getTumorRefCount(), variant.getTumorAltCount(), variant.getTumorReadDepth());

                    break;
                }
            }
        }
    }
}
