package com.hartwig.hmftools.isofox.neo;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.nextDown;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.utils.sv.SvRegion.positionsOverlap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.SINGLE_MAP_QUALITY;
import static com.hartwig.hmftools.isofox.common.GeneReadData.createGeneReadData;
import static com.hartwig.hmftools.isofox.common.ReadRecord.findOverlappingRegions;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.validExonMatch;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentMatcher.getNeoEpitopeSupport;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.neo.NeoEpitopeFile;
import com.hartwig.hmftools.common.utils.sv.SvRegion;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.BamSlicer;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.fusion.ReadGroup;

import org.jetbrains.annotations.NotNull;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;

public class NeoEpitopeReader
{
    private final IsofoxConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private final List<NeoEpitopeData> mNeoEpitopes;

    private final SamReader mSamReader;
    private final BamSlicer mBamSlicer;

    private GeneCollection mCurrentGenes;
    private NeoEpitopeData mCurrentNeoData;
    private final Map<String,ReadGroup> mReadGroups;

    public NeoEpitopeReader(final IsofoxConfig config, final EnsemblDataCache geneTransCache)
    {
        mConfig = config;
        mGeneTransCache = geneTransCache;

        mNeoEpitopes = Lists.newArrayList();
        mReadGroups = Maps.newHashMap();

        mCurrentGenes = null;
        mCurrentNeoData = null;

        mSamReader = mConfig.BamFile != null ?
                SamReaderFactory.makeDefault().referenceSequence(mConfig.RefGenomeFile).open(new File(mConfig.BamFile)) : null;

        boolean keepDuplicates = true;
        boolean keepSupplementaries = true;
        boolean keepSecondaries = true;
        int minMapQuality = keepSecondaries ? 0 : SINGLE_MAP_QUALITY;

        mBamSlicer = new BamSlicer(minMapQuality, keepDuplicates, keepSupplementaries, keepSecondaries);

        loadNeoEpitopes(mConfig.NeoEpitopeFile);
    }

    private void loadNeoEpitopes(final String filename)
    {
        if(filename == null || filename.isEmpty())
            return;

        try
        {
            final List<NeoEpitopeFile> sourceNEs = NeoEpitopeFile.read(filename);
            sourceNEs.forEach(x -> mNeoEpitopes.add(new NeoEpitopeData(x)));
            ISF_LOGGER.info("loaded {} neo-epitopes from file: {}", mNeoEpitopes.size(), filename);
        }
        catch(IOException exception)
        {
            ISF_LOGGER.error("failed to read neo-epitope file({})", filename, exception.toString());
        }
    }

    public void clearCache()
    {
        mReadGroups.clear();
    }

    public void calcFragmentSupport()
    {
        for(final NeoEpitopeData neData : mNeoEpitopes)
        {
            mCurrentNeoData = neData;
            clearCache();

            if(neData.isFusion())
            {
                calcFusionSupport(neData);
            }
            else
            {
                neData.setOrientation(mCurrentGenes.genes().get(0).GeneData.Strand);
                calcPointMutationSupport();
            }
        }

    }

    private void calcFusionSupport(final NeoEpitopeData neData)
    {
        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            initialiseGeneData(mCurrentNeoData.Source.GeneIds[fs]);

            int lowerPosition = min(mCurrentNeoData.Source.CodingBasePositions[fs], mCurrentNeoData.Positions[fs]);
            int upperPosition = max(mCurrentNeoData.Source.CodingBasePositions[fs], mCurrentNeoData.Positions[fs]);

            final SvRegion readRegion = new SvRegion(mCurrentNeoData.Chromosomes[fs], lowerPosition, upperPosition);

            mBamSlicer.slice(mSamReader, Lists.newArrayList(readRegion), this::processSamRecord);
        }

    }

    private void calcPointMutationSupport()
    {
        initialiseGeneData(mCurrentNeoData.Source.GeneIds[FS_UP]);

        final SvRegion readRegion = new SvRegion(
                mCurrentNeoData.Chromosomes[FS_UP],
                mCurrentNeoData.Source.CodingBasePositions[FS_UP], mCurrentNeoData.Source.CodingBasePositions[FS_DOWN]);


        mBamSlicer.slice(mSamReader, Lists.newArrayList(readRegion), this::processSamRecord);
    }

    private void initialiseGeneData(final String geneId)
    {
        final EnsemblGeneData geneData = mGeneTransCache.getGeneDataById(geneId);

        if(geneData == null)
        {
            ISF_LOGGER.error("gene({}) not found", geneId);
            return;
        }

        final List<GeneReadData> geneReadDataList = createGeneReadData(Lists.newArrayList(geneData), mGeneTransCache);
        mCurrentGenes = new GeneCollection(0, Lists.newArrayList(geneReadDataList));
    }

    private void processSamRecord(@NotNull final SAMRecord record)
    {
        final ReadRecord read = ReadRecord.from(record);

        read.processOverlappingRegions(findOverlappingRegions(mCurrentGenes.getExonRegions(), read));
        mCurrentGenes.setReadGeneCollections(read, mCurrentGenes.regionBounds());

        // only handle complete groups
        ReadGroup readGroup = mReadGroups.get(read.Id);

        if(readGroup == null)
        {
            mReadGroups.put(read.Id, new ReadGroup(read));
            return;
        }

        if(readGroup.isComplete())
        {
            processFragmentReads(readGroup);
            mReadGroups.remove(read.Id);
        }
        else
        {
            readGroup.Reads.add(read);
        }
    }

    private boolean isCandidateGroup(final ReadGroup readGroup)
    {
        if(!readGroup.isComplete())
            return false;

        // ignore unspliced reads
        if(mCurrentNeoData.isFusion())
        {
            boolean hasSplitRead = readGroup.Reads.stream().anyMatch(x -> x.containsSplit());

            if(!hasSplitRead && !readGroup.hasSuppAlignment())
                return false;
        }

        // check for support of the required transcript

        return true;
    }

    private void processFragmentReads(final ReadGroup readGroup)
    {
        if(!isCandidateGroup(readGroup))
            return;

        // each fragment must support at least one of the up & down transcripts by being fully exonic or matching an exon boundary
        // unspliced reads are skipped (in any of the transcripts in question)
        // for fusions, determine the gene covered by each read, and therefore the stream-edness of the read

        boolean[] hasSupportedTrans = { false, false };

        boolean singleGene = mCurrentNeoData.singleGene();

        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            if(singleGene && fs == FS_DOWN)
            {
                hasSupportedTrans[fs] = hasSupportedTrans[FS_UP];
                break;
            }

            for(String transName : mCurrentNeoData.Transcripts[fs])
            {
                boolean hasTrans = false;
                boolean hasSupport = true;

                for(ReadRecord read : readGroup.Reads)
                {
                    // any matched transcript needs to be fully exonic in every region, not just one
                    for(Map.Entry<RegionReadData,RegionMatchType> entry : read.getMappedRegions().entrySet())
                    {
                        final RegionReadData region = entry.getKey();

                        if(region.getTransExonRefs().stream().noneMatch(x -> transName.equals(x.TransName)))
                            continue;

                        hasTrans = true;

                        if(!validExonMatch(entry.getValue()))
                            hasSupport = false;
                    }
                }

                if(hasTrans && hasSupport)
                {
                    hasSupportedTrans[fs] = true;
                    break;
                }
            }
        }

        if(!hasSupportedTrans[FS_UP] || !hasSupportedTrans[FS_DOWN])
            return;

        // work out which of the sections of the neo-epitope may be supported
        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            if(singleGene && fs == FS_DOWN)
                break;

            final int[] codingBaseRange = mCurrentNeoData.getCodingBaseRange(fs);

            for(ReadRecord read : readGroup.Reads)
            {
                if(!read.Chromosome.equals(mCurrentNeoData.Chromosomes[fs]))
                    continue;

                // check that this read relates to the neo section
                if(!positionsOverlap(codingBaseRange[SE_START], codingBaseRange[SE_END], read.PosStart, read.PosEnd))
                    continue;

                NeoFragmentSupport support = getNeoEpitopeSupport(mCurrentNeoData, fs, read);
            }
        }
    }
}
