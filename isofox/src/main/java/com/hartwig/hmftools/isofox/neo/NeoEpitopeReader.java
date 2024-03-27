package com.hartwig.hmftools.isofox.neo;

import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_DOWN;
import static com.hartwig.hmftools.common.fusion.FusionCommon.FS_UP;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.IsofoxConstants.SINGLE_MAP_QUALITY;
import static com.hartwig.hmftools.isofox.common.GeneReadData.createGeneReadData;
import static com.hartwig.hmftools.isofox.common.ReadRecord.findOverlappingRegions;
import static com.hartwig.hmftools.isofox.common.RegionMatchType.validExonMatch;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentMatcher.checkBaseCoverage;
import static com.hartwig.hmftools.isofox.neo.NeoFragmentMatcher.findFusionSupport;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.neo.NeoEpitopeFile;
import com.hartwig.hmftools.common.neo.RnaNeoEpitope;
import com.hartwig.hmftools.common.bam.BamSlicer;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.isofox.IsofoxConfig;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.ReadRecord;
import com.hartwig.hmftools.isofox.common.RegionMatchType;
import com.hartwig.hmftools.isofox.common.RegionReadData;
import com.hartwig.hmftools.isofox.fusion.ChimericReadGroup;

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
    private final Map<String, ChimericReadGroup> mReadGroups;

    private BufferedWriter mWriter;

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

        boolean keepDuplicates = true; // will be dropped but need them to identify groups since supplementaries aren't marked
        boolean keepSupplementaries = true;
        boolean keepSecondaries = false;
        int minMapQuality = keepSecondaries ? 0 : SINGLE_MAP_QUALITY;

        mBamSlicer = new BamSlicer(minMapQuality, keepDuplicates, keepSupplementaries, keepSecondaries);

        String neoepitopeFile = NeoEpitopeFile.generateFilename(mConfig.NeoDir, mConfig.SampleId);
        loadNeoEpitopes(neoepitopeFile);
        initialiseWriter();
    }

    private void loadNeoEpitopes(final String filename)
    {
        if(filename == null || filename.isEmpty())
            return;

        try
        {
            final List<NeoEpitopeFile> sourceNEs = NeoEpitopeFile.read(filename);
            sourceNEs.forEach(x -> mNeoEpitopes.add(new NeoEpitopeData(x)));
            ISF_LOGGER.info("sample({}) loaded {} neo-epitopes from file: {}", mConfig.SampleId, mNeoEpitopes.size(), filename);
        }
        catch(IOException exception)
        {
            ISF_LOGGER.error("failed to read neo-epitope file({}): {}", filename, exception.toString());
        }
    }

    public void clearCache()
    {
        mReadGroups.clear();
    }

    private boolean filterOnRestrictedGenes(final NeoEpitopeData neData)
    {
        if(mConfig.Filters.RestrictedGeneIds.isEmpty())
            return true;

        return mConfig.Filters.RestrictedGeneIds.contains(neData.Source.GeneIds[FS_UP])
                && mConfig.Filters.RestrictedGeneIds.contains(neData.Source.GeneIds[FS_DOWN]);
    }

    public void calcFragmentSupport()
    {
        for(final NeoEpitopeData neData : mNeoEpitopes)
        {
            clearCache();

            if(!filterOnRestrictedGenes(neData))
                continue;

            mCurrentNeoData = neData;

            if(neData.isPointMutation())
                continue;

            calcFusionSupport();

            writeData(neData);
        }

        closeBufferedWriter(mWriter);
    }

    private void calcFusionSupport()
    {
        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            initialiseGeneData(mCurrentNeoData.Source.GeneIds[fs]);

            final ChrBaseRegion readRegion = new ChrBaseRegion(mCurrentNeoData.Chromosomes[fs], mCurrentNeoData.Source.CodingBasePositions[fs]);
            mBamSlicer.slice(mSamReader, readRegion, this::processSamRecord);
        }

        mReadGroups.values().forEach(x -> processFragmentReads(x));
    }

    private void initialiseGeneData(final String geneId)
    {
        mCurrentGenes = null;

        final GeneData geneData = mGeneTransCache.getGeneDataById(geneId);

        if(geneData == null)
        {
            ISF_LOGGER.error("gene({}) not found", geneId);
            return;
        }

        final List<GeneReadData> geneReadDataList = createGeneReadData(Lists.newArrayList(geneData), mGeneTransCache);
        mCurrentGenes = new GeneCollection(0, Lists.newArrayList(geneReadDataList));

        if(mConfig.RefGenomeFile != null)
        {
            for(RegionReadData region : mCurrentGenes.getExonRegions())
            {
                final String regionRefBases = mConfig.RefGenome.getBaseString(region.chromosome(), region.start(), region.end());
                region.setRefBases(regionRefBases);
            }
        }
    }

    private void processSamRecord(final SAMRecord record)
    {
        final ReadRecord read = ReadRecord.from(record);

        read.processOverlappingRegions(findOverlappingRegions(mCurrentGenes.getExonRegions(), read));
        mCurrentGenes.setReadGeneCollections(read, mCurrentGenes.regionBounds());

        // only handle complete groups
        ChimericReadGroup readGroup = mReadGroups.get(read.Id);

        if(readGroup == null)
        {
            mReadGroups.put(read.Id, new ChimericReadGroup(read));
            return;
        }

        readGroup.addRead(read);

        if(readGroup.isComplete())
        {
            processFragmentReads(readGroup);
            mReadGroups.remove(read.Id);
        }
    }

    private boolean isCandidateGroup(final ChimericReadGroup readGroup)
    {
        // ignore unspliced reads
        if(mCurrentNeoData.isFusion())
        {
            if(readGroup.hasSuppAlignment())
                return true;

            return readGroup.reads().stream().anyMatch(x -> x.containsSplit());
        }

        // check for support of the required transcript

        return true;
    }

    private void processFragmentReads(final ChimericReadGroup readGroup)
    {
        if(readGroup.reads().stream().anyMatch(x -> x.isDuplicate()))
            return;

        checkBaseCoverage(mCurrentNeoData, readGroup);

        if(!isCandidateGroup(readGroup))
            return;

        // each fragment must support at least one of the up & down transcripts by being fully exonic or matching an exon boundary
        // unspliced reads are skipped (in any of the transcripts in question)
        // for fusions, determine the gene covered by each read, and therefore the stream-edness of the read

        boolean[] hasSupportedTrans = { false, false };

        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            if(mCurrentNeoData.singleGene() && fs == FS_DOWN)
            {
                hasSupportedTrans[fs] = hasSupportedTrans[FS_UP];
                break;
            }

            for(String transName : mCurrentNeoData.Transcripts[fs])
            {
                boolean hasTrans = false;
                boolean hasSupport = true;

                for(ReadRecord read : readGroup.reads())
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
        NeoFragmentSupport readGroupSupport = new NeoFragmentSupport();

        for(int fs = FS_UP; fs <= FS_DOWN; ++fs)
        {
            final int[] codingBaseRange = mCurrentNeoData.getCodingBaseRange(fs);

            for(ReadRecord read : readGroup.reads())
            {
                if(!read.Chromosome.equals(mCurrentNeoData.Chromosomes[fs]))
                    continue;

                // check that this read covers some part of the neo section
                if(read.getMappedRegionCoords(false).stream()
                        .noneMatch(x -> positionsOverlap(codingBaseRange[SE_START], codingBaseRange[SE_END], x[SE_START], x[SE_END])))
                {
                    continue;
                }

                NeoFragmentSupport support = findFusionSupport(mCurrentNeoData, fs, read);
                readGroupSupport.setMax(support);

                ISF_LOGGER.trace("neo({}:{}) support({}) read({})",
                        mCurrentNeoData.Source.Id, mCurrentNeoData.Source.VariantInfo, support, readGroup.id());

                /*
                if(support.NovelFragments[EXACT_MATCH] > 0)
                {
                    checkBaseCoverage(mCurrentNeoData, readGroup);
                }
                */
            }
        }

        mCurrentNeoData.getFragmentSupport().combineSupport(readGroupSupport);
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = RnaNeoEpitope.generateFilename(mConfig.OutputDir, mConfig.SampleId);

            mWriter = createBufferedWriter(outputFileName, false);
            mWriter.write(RnaNeoEpitope.header());
            mWriter.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to create neo-epitope writer: {}", e.toString());
        }
    }

    private void writeData(final NeoEpitopeData neData)
    {
        try
        {
            mWriter.write(RnaNeoEpitope.toString(neData.asRnaFile()));
            mWriter.newLine();
        }
        catch (IOException e)
        {
            ISF_LOGGER.error("failed to write neo-epitope data: {}", e.toString());
        }
    }
}
