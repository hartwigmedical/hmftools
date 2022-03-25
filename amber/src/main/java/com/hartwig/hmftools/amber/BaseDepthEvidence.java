package com.hartwig.hmftools.amber;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.amber.AmberSite;
import com.hartwig.hmftools.common.amber.BaseDepth;
import com.hartwig.hmftools.common.amber.BaseDepthFactory;
import com.hartwig.hmftools.common.amber.ModifiableBaseDepth;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelector;
import com.hartwig.hmftools.common.genome.position.GenomePositionSelectorFactory;
import com.hartwig.hmftools.common.genome.position.GenomePositions;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.samtools.BamSlicer;
import com.hartwig.hmftools.common.utils.collection.Multimaps;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class BaseDepthEvidence implements Callable<BaseDepthEvidence>
{
    private final String mChromosome;
    private final String mBamFile;
    private final SamReaderFactory mSamReaderFactory;
    private final List<ModifiableBaseDepth> mEvidenceMap;
    private final GenomePositionSelector<ModifiableBaseDepth> mSelector;
    private final BaseDepthFactory mBafFactory;
    private final BamSlicer mBamSlicer;
    private final List<GenomeRegion> mBafRegions;

    public BaseDepthEvidence(
            int typicalReadDepth, int minMappingQuality, int minBaseQuality, final String contig, final String bamFile,
            final SamReaderFactory samReaderFactory, final List<AmberSite> bafLoci)
    {
        mBafFactory = new BaseDepthFactory(minBaseQuality);
        mChromosome = contig;
        mBamFile = bamFile;
        mSamReaderFactory = samReaderFactory;

        mBafRegions = GenomeRegions.fromSortedGenomePositions(typicalReadDepth, bafLoci);

        AMB_LOGGER.debug("created {} regions from {} baf loci", mBafRegions.size(), bafLoci.size());

        mEvidenceMap = bafLoci.stream().map(BaseDepthFactory::fromAmberSite).collect(Collectors.toList());
        mSelector = GenomePositionSelectorFactory.create(Multimaps.fromPositions(mEvidenceMap));
        mBamSlicer = new BamSlicer(minMappingQuality);
    }

    @NotNull
    public String contig()
    {
        return mChromosome;
    }

    @NotNull
    public List<BaseDepth> evidence()
    {
        return new ArrayList<>(mEvidenceMap);
    }

    @Override
    public BaseDepthEvidence call() throws Exception
    {
        AMB_LOGGER.debug("germline bam processing {} regions", mBafRegions.size());
        try(SamReader reader = mSamReaderFactory.open(new File(mBamFile)))
        {
            mBamSlicer.sliceNoDups(reader, mBafRegions, this::record);
        }
        AMB_LOGGER.debug("germline bam processing done");

        return this;
    }

    private void record(@NotNull final SAMRecord record)
    {
        mSelector.select(asRegion(record), bafEvidence -> mBafFactory.addEvidence(bafEvidence, record));
    }

    @NotNull
    private static GenomeRegion asRegion(@NotNull final SAMRecord record)
    {
        return GenomeRegions.create(record.getContig(), record.getAlignmentStart(), record.getAlignmentEnd());
    }

    private void processRecord(@NotNull final SAMRecord record)
    {
        if (!record.getContig().equals(mChromosome))
            return;

        // use binary search on the evidence list
        // we rely on the fact that they all come from same chromosome
        int lowerBound = Collections.binarySearch(mEvidenceMap, GenomePositions.create(record.getContig(), record.getAlignmentStart()),
                Comparator.comparingInt(GenomePosition::position));

        if (lowerBound < 0)
        {
            // if can't find it then binarySearch returns -(insertion point) - 1
            // so x = -(insertion point) - 1
            //  insertion point = -x - 1
            lowerBound = -lowerBound - 1;
        }

        for (int i = lowerBound; i < mEvidenceMap.size(); ++i)
        {
            ModifiableBaseDepth e = mEvidenceMap.get(i);

            if (e.position() < record.getAlignmentStart())
            {
                continue;
            }
            if (e.position() > record.getAlignmentEnd())
            {
                break;
            }

            // here includes all positions in [align start, align end]
            mBafFactory.addEvidence(e, record);
        }
    }


    // compare genome regions
    void compareRegions(final List<GenomeRegion> gpList1, final List<GenomeRegion> gpList2)
    {
        var itr1 = gpList1.iterator();
        var itr2 = gpList2.iterator();

        while (itr1.hasNext() && itr2.hasNext())
        {
            var gp1 = itr1.next();
            var gp2 = itr2.next();

            if (!gp1.chromosome().equals(gp2.chromosome()) ||
                gp1.start() != gp2.start() ||
                gp1.end() != gp2.end())
            {
                AMB_LOGGER.info("regions mismatch!!");
            }
        }

        if (itr1.hasNext())
        {
            var gp1 = itr1.next();
            AMB_LOGGER.info("regions mismatch!!");
        }
        if  (itr2.hasNext())
        {
            var gp2 = itr2.next();
            AMB_LOGGER.info("regions mismatch!!");
        }

        AMB_LOGGER.info("regions match!!");
    }
}
