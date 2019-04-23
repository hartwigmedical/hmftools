package com.hartwig.hmftools.bachelor;

import static com.hartwig.hmftools.bachelor.BachelorPostProcess.REF_GENOME;

import java.io.File;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.bachelor.types.BachelorGermlineVariant;
import com.hartwig.hmftools.common.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.hotspot.VariantHotspot;
import com.hartwig.hmftools.common.hotspot.VariantHotspotEvidence;
import com.hartwig.hmftools.common.hotspot.VariantHotspotEvidenceFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;

public class BamCountReader
{
    private IndexedFastaSequenceFile mIndexedFastaSequenceFile;
    SamReader mTumorReader;

    private static String TUMOR_BAM_PATH = "tumor_bam_path";
    private static int DEFAULT_MIN_BASE_QUALITY = 13;
    private static int DEFAULT_MIN_MAPPING_QUALITY = 1;

    private static final Logger LOGGER = LogManager.getLogger(BamCountReader.class);

    public BamCountReader()
    {
        mIndexedFastaSequenceFile = null;
        mTumorReader = null;
    }

    public static void addCmdLineOptions(Options options)
    {
        options.addOption(TUMOR_BAM_PATH, true, "Max batch directories to batch process");
    }

    public boolean initialise(final CommandLine cmd, IndexedFastaSequenceFile indexedFastaSequenceFile)
    {
        mIndexedFastaSequenceFile = indexedFastaSequenceFile;
        String bamPath = cmd.getOptionValue(TUMOR_BAM_PATH);
        File refGenomeFile = new File(cmd.getOptionValue(REF_GENOME));
        mTumorReader = SamReaderFactory.makeDefault().referenceSequence(refGenomeFile).open(new File(bamPath));
        return true;
    }

    public void readBamCounts(List<BachelorGermlineVariant> bachRecords)
    {
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

        if(tumorEvidence.size() != bachRecords.size())
        {
            LOGGER.error("incomplete BAM evidence read");
            return;
        }

        for(int i = 0; i < bachRecords.size(); ++i)
        {
            BachelorGermlineVariant variant = bachRecords.get(i);
            VariantHotspotEvidence evidence = tumorEvidence.get(i);

            variant.setTumorData(evidence.altSupport(), evidence.readDepth());

            LOGGER.debug("chr({}) position({}) matched, counts(ref={} alt={} depth={})",
                    variant.Chromosome, variant.Position,
                    variant.getTumorRefCount(), variant.getTumorAltCount(), variant.getTumorReadDepth());

        }
    }

}
