package com.hartwig.hmftools.isofox.novel.cohort;

import static java.lang.Math.abs;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionContext.SPLICE_JUNC;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.EXON_INTRON;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.MIXED_TRANS;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.NOVEL_3_PRIME;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.NOVEL_5_PRIME;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.NOVEL_EXON;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.SKIPPED_EXONS;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionType.UNKNOWN;
import static com.hartwig.hmftools.common.utils.Strings.appendStrList;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionsWithin;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.PASS_FILTER;
import static com.hartwig.hmftools.isofox.IsofoxConfig.ISF_LOGGER;
import static com.hartwig.hmftools.isofox.cohort.AnalysisType.ALT_SPLICE_JUNCTION;
import static com.hartwig.hmftools.isofox.cohort.CohortConfig.formSampleFilenames;
import static com.hartwig.hmftools.isofox.novel.cohort.AcceptorDonorType.ACCEPTOR;
import static com.hartwig.hmftools.isofox.novel.cohort.AcceptorDonorType.DONOR;
import static com.hartwig.hmftools.isofox.novel.cohort.AcceptorDonorType.NONE;
import static com.hartwig.hmftools.isofox.novel.cohort.AltSjCohortAnalyser.loadFile;
import static com.hartwig.hmftools.isofox.novel.cohort.SpliceSiteCache.PSI_NO_RATE;
import static com.hartwig.hmftools.isofox.novel.cohort.SpliceSiteCache.SS_SUPPORT;
import static com.hartwig.hmftools.isofox.novel.cohort.SpliceSiteCache.SS_TRAVERSED;
import static com.hartwig.hmftools.isofox.novel.cohort.SpliceSiteCache.calcSupportRate;
import static com.hartwig.hmftools.isofox.novel.cohort.SpliceVariantMatchType.DISRUPTION;
import static com.hartwig.hmftools.isofox.novel.cohort.SpliceVariantMatchType.NOVEL;
import static com.hartwig.hmftools.isofox.results.ResultsWriter.ITEM_DELIM;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.ExonData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionContext;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionType;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.sv.StructuralVariantData;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.isofox.cohort.CohortConfig;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Options;
import org.jooq.Record;
import org.jooq.Result;

public class SpliceVariantMatcher
{
    private final CohortConfig mConfig;
    private final EnsemblDataCache mGeneTransCache;

    private BufferedWriter mWriter;
    private final Map<String,Integer> mFieldsMap;

    private final Map<String,EnsemblGeneData> mGeneDataMap;
    private final AltSjFilter mAltSjFilter;

    private final SpliceVariantCache mDataCache;
    private final SpliceSiteCache mSpliceSiteCache;
    private final List<SpliceVariantMatchType> mMatchTypes;

    protected static final String SOMATIC_VARIANT_FILE = "somatic_variant_file";
    protected static final String SV_BREAKEND_FILE = "sv_breakend_file";
    protected static final String COHORT_ALT_SJ_FILE = "cohort_alt_sj_file";

    private static final String INCLUDE_ALL_TRANSCRIPTS = "include_all_transcripts";
    private static final String SPLICE_VARIANT_TYPES = "splice_variant_types";
    protected static final String WRITE_VARIANT_CACHE = "write_variant_cache";

    private static final int MIN_ALT_SJ_LENGTH = 65;
    private static final int SPLICE_REGION_NON_CODING_DISTANCE = 10;
    private static final int SPLICE_REGION_CODING_DISTANCE = 2;
    private static final int CLOSE_ALT_SJ_DISTANCE = 5;
    private static final int ALT_SJ_MISALIGN_BUFFER = 5;

    public SpliceVariantMatcher(final CohortConfig config, final CommandLine cmd)
    {
        mConfig = config;
        mFieldsMap = Maps.newHashMap();
        mGeneDataMap = Maps.newHashMap();

        boolean allTranscripts = cmd.hasOption(INCLUDE_ALL_TRANSCRIPTS);

        mGeneTransCache = new EnsemblDataCache(mConfig.EnsemblDataCache, RefGenomeVersion.V37);
        mGeneTransCache.setRequiredData(true, false, false, !allTranscripts);
        mGeneTransCache.load(false);

        mAltSjFilter = new AltSjFilter(mConfig.RestrictedGeneIds, mConfig.ExcludedGeneIds, 0);
        mWriter = null;

        mMatchTypes = Lists.newArrayList();

        if(cmd.hasOption(SPLICE_VARIANT_TYPES))
        {
            Arrays.stream(cmd.getOptionValue(SPLICE_VARIANT_TYPES).split(";", -1))
                    .map(x -> SpliceVariantMatchType.valueOf(x)).forEach(x -> mMatchTypes.add(x));
        }

        initialiseWriter();

        mDataCache = new SpliceVariantCache(config, cmd);

        for(String geneId : mConfig.RestrictedGeneIds)
        {
            final EnsemblGeneData geneData = mGeneTransCache.getGeneDataById(geneId);

            if(geneData != null)
            {
                mGeneDataMap.put(geneData.GeneName, geneData);
            }
        }

        mSpliceSiteCache = new SpliceSiteCache(config, cmd);
    }

    public static void addCmdLineOptions(final Options options)
    {
        options.addOption(SOMATIC_VARIANT_FILE, true, "File with somatic variants potentially affecting splicing");
        options.addOption(SV_BREAKEND_FILE, true, "File with cached SV positions");
        options.addOption(COHORT_ALT_SJ_FILE, true, "Cohort frequency for alt SJs");
        options.addOption(WRITE_VARIANT_CACHE, false, "Write out somatic variants for subsequent non-DB loading");
        options.addOption(INCLUDE_ALL_TRANSCRIPTS, false, "Consider all transcripts, not just canonical (default: false)");
        options.addOption(SPLICE_VARIANT_TYPES, true, "Types to include in analysis: NOVEL, DISRUPTION, empty=ALL");
    }

    public void processAltSpliceJunctions()
    {
        if(mGeneTransCache == null)
        {
            ISF_LOGGER.error("missing Ensembl data cache,");
            return;
        }

        if(!mDataCache.hasCachedSomaticVariants() && mConfig.DbAccess == null)
        {
            ISF_LOGGER.error("missing DB connection or cached somatic variants,");
            return;
        }

        final List<Path> filenames = Lists.newArrayList();

        if(!formSampleFilenames(mConfig, ALT_SPLICE_JUNCTION, filenames))
            return;

        for(int i = 0; i < mConfig.SampleData.SampleIds.size(); ++i)
        {
            final String sampleId = mConfig.SampleData.SampleIds.get(i);
            final Path altSJFile = filenames.get(i);

            final List<AltSpliceJunctionFile> altSJs = loadFile(altSJFile, mFieldsMap, mAltSjFilter);

            ISF_LOGGER.debug("{}: sample({}) loaded {} alt-SJ records", i, sampleId, altSJs.size());
            evaluateSpliceVariants(sampleId, altSJs);
        }

        ISF_LOGGER.info("splice variant matching complete");

        closeBufferedWriter(mWriter);mDataCache.close();
    }

    private void evaluateSpliceVariants(final String sampleId, final List<AltSpliceJunctionFile> altSpliceJunctions)
    {
        final List<SpliceVariant> spliceVariants = getSomaticVariants(sampleId);

        if(spliceVariants == null || spliceVariants.isEmpty())
        {
            ISF_LOGGER.debug("sampleId({}) no splice variants found", sampleId);
            return;
        }

        final Map<String,List<Integer>> svBreakends = getStructuralVariants(sampleId);

        final Map<String,List<AltSpliceJunctionFile>> candidateAltSJs = Maps.newHashMap();

        for(AltSpliceJunctionFile altSJ : altSpliceJunctions)
        {
            if(altSJ.length() < MIN_ALT_SJ_LENGTH)
                continue;

            if(mMatchTypes.contains(NOVEL))
            {
                if(!validCrypticType(altSJ.Type))
                    continue;
            }

            List<Integer> breakends = svBreakends.get(altSJ.Chromosome);

            if(breakends != null)
            {
                if(breakends.stream().anyMatch(x -> positionWithin(x, altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END])))
                    continue;
            }

            List<AltSpliceJunctionFile> chrAltSJs = candidateAltSJs.get(altSJ.Chromosome);

            if(chrAltSJs == null)
            {
                chrAltSJs = Lists.newArrayList(altSJ);
                candidateAltSJs.put(altSJ.Chromosome, chrAltSJs);
            }
            else
            {
                chrAltSJs.add(altSJ);
            }
        }

        ISF_LOGGER.debug("sampleId({}) evaluating {} splice variants vs altSJs({})",
                sampleId, spliceVariants.size(), candidateAltSJs.values().stream().mapToInt(x -> x.size()).sum());

        mSpliceSiteCache.loadSampleSpliceSites(sampleId);

        spliceVariants.forEach(x -> evaluateSpliceVariant(sampleId, x, candidateAltSJs));

        mDataCache.writeVariantCache(sampleId, spliceVariants, svBreakends);
    }

    private List<SpliceVariant> getSomaticVariants(final String sampleId)
    {
        if(mDataCache.hasCachedSomaticVariants())
            return mDataCache.retrieveSomaticVariants(sampleId);

        final List<SpliceVariant> spliceVariants = Lists.newArrayList();

        final Result<Record> result = mConfig.DbAccess.context().select().from(Tables.SOMATICVARIANT)
                .where(Tables.SOMATICVARIANT.SAMPLEID.eq(sampleId))
                .and(Tables.SOMATICVARIANT.FILTER.eq(PASS_FILTER))
                .and(Tables.SOMATICVARIANT.GENE.notEqual(""))
                .fetch();

        for (Record record : result)
        {
            // filter to specific gene list
            final String gene = record.getValue(Tables.SOMATICVARIANT.GENE);

            if(!mGeneDataMap.isEmpty() && !mGeneDataMap.containsKey(gene))
                continue;

            String chromosome = record.getValue(Tables.SOMATICVARIANT.CHROMOSOME);
            int position = record.getValue(Tables.SOMATICVARIANT.POSITION);
            String ref = record.getValue(Tables.SOMATICVARIANT.REF);
            String alt = record.getValue(Tables.SOMATICVARIANT.ALT);
            Integer localPhaseSet = record.get(Tables.SOMATICVARIANT.LOCALPHASESET);
            String triNecContext = record.getValue(SOMATICVARIANT.TRINUCLEOTIDECONTEXT);
            String canonicalEffect = record.getValue(SOMATICVARIANT.CANONICALEFFECT);
            String canonicalCodingImpact = record.getValue(SOMATICVARIANT.CANONICALHGVSCODINGIMPACT);
            VariantType type = VariantType.valueOf(record.getValue(SOMATICVARIANT.TYPE));

            spliceVariants.add(new SpliceVariant(
                    gene, chromosome, position, type, ref, alt,
                    canonicalEffect, canonicalCodingImpact, triNecContext, localPhaseSet != null ? localPhaseSet : -1));
        }

        return spliceVariants;
    }

    private Map<String,List<Integer>> getStructuralVariants(final String sampleId)
    {
        if(mDataCache.hasCachedSvBreakends())
            return mDataCache.retrieveSvBreakends(sampleId);

        final Map<String,List<Integer>> svBreakends = Maps.newHashMap();

        if(mConfig.DbAccess == null)
            return svBreakends;

        final List<StructuralVariantData> structuralVariants = mConfig.DbAccess.readStructuralVariantData(sampleId);

        for(StructuralVariantData sv : structuralVariants)
        {
            List<Integer> positions = svBreakends.get(sv.startChromosome());

            if(positions == null)
                svBreakends.put(sv.startChromosome(), Lists.newArrayList(sv.startPosition()));
            else
                positions.add(sv.startPosition());

            if(sv.type() != StructuralVariantType.INF && sv.type() != StructuralVariantType.SGL)
            {
                positions = svBreakends.get(sv.endChromosome());

                if(positions == null)
                    svBreakends.put(sv.endChromosome(), Lists.newArrayList(sv.endPosition()));
                else
                    positions.add(sv.endPosition());
            }
        }

        return svBreakends;
    }

    private void evaluateSpliceVariant(
            final String sampleId, final SpliceVariant variant, final Map<String,List<AltSpliceJunctionFile>> altSpliceJunctions)
    {
        final EnsemblGeneData geneData = mGeneTransCache.getGeneDataByName(variant.GeneName);

        if(geneData == null)
        {
            ISF_LOGGER.error("variant({}) gene data not found", variant);
            return;
        }

        List<AltSpliceJunctionFile> chrAltSJs = altSpliceJunctions.get(variant.Chromosome);

        if(chrAltSJs == null || chrAltSJs.isEmpty())
            return;

        final List<AltSpliceJunctionFile> matchedAltSJs = Lists.newArrayList();

        if(mMatchTypes.isEmpty() || mMatchTypes.contains(NOVEL))
            matchedAltSJs.addAll(findCrypticAltSpliceJunctions(sampleId, chrAltSJs, variant, geneData));

        if(mMatchTypes.isEmpty() || mMatchTypes.contains(DISRUPTION))
            findRelatedAltSpliceJunctions(sampleId, chrAltSJs, variant, geneData, matchedAltSJs);
    }

    private static boolean withinRange(int pos1, int pos2, int distance) { return abs(pos1 - pos2) <= distance; }

    private static boolean validDisruptionType(AltSpliceJunctionType type, AcceptorDonorType acceptorDonorType)
    {
        if(type == NOVEL_5_PRIME)
            return acceptorDonorType == DONOR;

        if(type == NOVEL_3_PRIME)
            return acceptorDonorType == ACCEPTOR;

        return type == SKIPPED_EXONS || type == NOVEL_EXON || type == MIXED_TRANS;
    }

    private void findRelatedAltSpliceJunctions(
            final String sampleId, final List<AltSpliceJunctionFile> altSpliceJunctions,
            final SpliceVariant variant, final EnsemblGeneData geneData, final List<AltSpliceJunctionFile> matchedAltSJs)
    {
        final List<TranscriptData> transDataList = mGeneTransCache.getTranscripts(geneData.GeneId);

        boolean isWithinSpliceRegion = false;

        AcceptorDonorType closestAcceptorDonorType = AcceptorDonorType.NONE;
        int closestExonDistance = 0;
        int closestExonPosition = 0;
        int exonPosition = 0;
        String closestTransStr = "";

        for(final TranscriptData transData : transDataList)
        {
            for(int i = 0; i < transData.exons().size(); ++i)
            {
                final ExonData exon = transData.exons().get(i);

                AcceptorDonorType acceptorDonorType = AcceptorDonorType.NONE;
                int exonDistance = 0;

                if((i == 0 && withinRange(exon.Start, variant.Position, SPLICE_REGION_NON_CODING_DISTANCE))
                || (i == transData.exons().size() - 1 && withinRange(exon.End, variant.Position, SPLICE_REGION_NON_CODING_DISTANCE)))
                {
                    // ignore matching at the first or last base of a transcript
                    continue;
                }

                if(positionWithin(variant.Position, exon.Start, exon.End))
                {
                    // variant within an exon
                    if(withinRange(exon.Start, variant.Position, SPLICE_REGION_CODING_DISTANCE))
                    {
                        acceptorDonorType = geneData.Strand == POS_STRAND ? ACCEPTOR : DONOR;
                        exonDistance = -abs(exon.Start - variant.Position) - 1; // first base inside exon is -1 by convention
                        exonPosition = exon.Start;
                    }
                    else if(withinRange(exon.End, variant.Position, SPLICE_REGION_CODING_DISTANCE))
                    {
                        acceptorDonorType = geneData.Strand == POS_STRAND ? DONOR : ACCEPTOR;
                        exonDistance = -abs(exon.End - variant.Position) - 1;
                        exonPosition = exon.End;
                    }
                    else
                    {
                        continue;
                    }
                }
                else
                {
                    if(withinRange(exon.Start, variant.Position, SPLICE_REGION_NON_CODING_DISTANCE))
                    {
                        acceptorDonorType = geneData.Strand == POS_STRAND ? ACCEPTOR : DONOR;
                        exonDistance = abs(exon.Start - variant.Position);
                        exonPosition = exon.Start;
                    }
                    else if(withinRange(exon.End, variant.Position, SPLICE_REGION_NON_CODING_DISTANCE))
                    {
                        acceptorDonorType = geneData.Strand == POS_STRAND ? DONOR : ACCEPTOR;
                        exonDistance = abs(exon.End - variant.Position);
                        exonPosition = exon.End;
                    }
                    else
                    {
                        continue;
                    }
                }

                isWithinSpliceRegion = true;

                if(closestAcceptorDonorType == NONE || abs(exonDistance) < abs(closestExonDistance))
                {
                    closestAcceptorDonorType = acceptorDonorType;
                    closestExonDistance = exonDistance;
                    closestExonPosition = exonPosition;
                    closestTransStr = String.format("%s;%d", transData.TransName, exon.Rank);
                }

                final ExonData prevExon = i > 0 ? transData.exons().get(i - 1) : null;
                final ExonData nextExon = i < transData.exons().size() - 1 ? transData.exons().get(i + 1) : null;

                // look for any alt SJs which match with the somatic variant
                for(final AltSpliceJunctionFile altSJ : altSpliceJunctions)
                {
                    if(matchedAltSJs.contains(altSJ))
                        continue;

                    if(!altSJ.GeneId.equals(geneData.GeneId))
                        continue;

                    if(!validDisruptionType(altSJ.Type, acceptorDonorType))
                        continue;

                    if(altSJ.Type == SKIPPED_EXONS)
                    {
                        // SKIPPED-EXONS - must match the previous and next exons, ie only skips the exon in question
                        if(prevExon == null || nextExon == null)
                            continue;

                        if(altSJ.SpliceJunction[SE_START] != prevExon.End || altSJ.SpliceJunction[SE_END] != nextExon.Start)
                            continue;
                    }
                    else
                    {
                        int minPosition = prevExon != null ? prevExon.End : exon.Start;
                        int maxPosition = nextExon != null ? nextExon.Start : exon.End;

                        // allow a margin for misalignment
                        minPosition -= ALT_SJ_MISALIGN_BUFFER;
                        maxPosition += ALT_SJ_MISALIGN_BUFFER;

                        if(!positionsWithin(altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END], minPosition, maxPosition))
                            continue;
                    }

                    ISF_LOGGER.trace("sampleId({}) variant({}:{}) gene({}) transcript({}) matched altSJ({})",
                            sampleId, variant.Chromosome, variant.Position, geneData.GeneName, transData.TransName, altSJ.toString());

                    final String transDataStr = String.format("%s;%d", transData.TransName, exon.Rank);
                    writeMatchData(sampleId, variant, geneData, altSJ, DISRUPTION, acceptorDonorType, exonDistance, exonPosition, transDataStr);

                    matchedAltSJs.add(altSJ);
                }
            }
        }

        if(matchedAltSJs.isEmpty() && isWithinSpliceRegion)
        {
            writeMatchData(
                    sampleId, variant, geneData, null, DISRUPTION,
                    closestAcceptorDonorType, closestExonDistance, closestExonPosition, closestTransStr);
        }
    }

    private static boolean validCrypticType(AltSpliceJunctionType type)
    {
        return type == NOVEL_5_PRIME || type == NOVEL_3_PRIME || type == NOVEL_EXON || type == EXON_INTRON;
    }

    private List<AltSpliceJunctionFile> findCrypticAltSpliceJunctions(
            final String sampleId, final List<AltSpliceJunctionFile> altSpliceJunctions, final SpliceVariant variant,
            final EnsemblGeneData geneData)
    {
        final List<AltSpliceJunctionFile> closeAltSJs = Lists.newArrayList();

        for(AltSpliceJunctionFile altSJ : altSpliceJunctions)
        {
            if(!validCrypticType(altSJ.Type))
                continue;

            Integer closeIndex = null;

            for(int se = SE_START; se <= SE_END; ++se)
            {
                if(!withinRange(altSJ.SpliceJunction[se], variant.Position, CLOSE_ALT_SJ_DISTANCE))
                    continue;

                if(altSJ.RegionContexts[se] == SPLICE_JUNC)
                    continue;

                closeIndex = se;
                break;
            }

            if(closeIndex == null)
                continue;

            List<String> allTransNames = Lists.newArrayList();

            for(int se = SE_START; se <= SE_END; ++se)
            {
                final String[] transNames = altSJ.TranscriptNames[se].split(ITEM_DELIM);
                for(String transName : transNames)
                {
                    if(!transName.equals("NONE") && !allTransNames.contains(transName))
                        allTransNames.add(transName);
                }
            }

            // acceptor donor based on gene strand and which end of the splice junction is matched
            AcceptorDonorType acceptorDonorType = ((closeIndex == SE_START) == (geneData.Strand == POS_STRAND)) ? DONOR : ACCEPTOR;

            writeMatchData(
                    sampleId, variant, geneData, altSJ, NOVEL, acceptorDonorType, -1, null,
                    appendStrList(allTransNames, ITEM_DELIM.charAt(0)));

            closeAltSJs.add(altSJ);
        }

        return closeAltSJs;
    }

    private String getBaseContext(final SpliceVariant variant, int splicePosition, final AcceptorDonorType acceptorDonorType)
    {
        if(mConfig.RefGenome == null)
            return "";

        return SpliceVariant.getBaseContext(
                variant.Chromosome, variant.Position, variant.Ref, variant.Alt, splicePosition, acceptorDonorType, mConfig.RefGenome);
    }

    private void initialiseWriter()
    {
        try
        {
            final String outputFileName = mConfig.formCohortFilename("splice_variant_matching.csv");
            mWriter = createBufferedWriter(outputFileName, false);

            mWriter.write("SampleId,MatchType,GeneId,GeneName,Chromosome,Strand");
            mWriter.write(",Position,Type,CodingEffect,Ref,Alt,HgvsImpact,TriNucContext,LocalPhaseSet,AccDonType,ExonDistance,ExonPosition");
            mWriter.write(",AsjType,AsjContextStart,AsjContextEnd,AsjPosStart,AsjPosEnd,FragmentCount,DepthStart,DepthEnd");
            mWriter.write(",BasesStart,BasesEnd,AsjCohortCount,AsjTransData");
            mWriter.write(",PsiSample,PsiSampleTotalFrags,PsiCohortPerc,PsiCohortAvg,PsiCohortTotalFrags");

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write splice variant file: {}", e.toString());
        }
    }

    private void writeMatchData(
            final String sampleId, final SpliceVariant variant, final EnsemblGeneData geneData, final AltSpliceJunctionFile altSJ,
            final SpliceVariantMatchType matchType, final AcceptorDonorType accDonType,
            int exonBaseDistance, Integer exonPosition, final String transDataStr)
    {
        try
        {
            mWriter.write(String.format("%s,%s,%s,%s,%s,%d",
                    sampleId, matchType, geneData.GeneId, geneData.GeneName, geneData.Chromosome, geneData.Strand));

            mWriter.write(String.format(",%d,%s,%s,%s,%s,%s,%s,%d,%s,%d,%d",
                    variant.Position, variant.Type, variant.CodingEffect, variant.Ref, variant.Alt,
                    variant.CodingImpact, variant.TriNucContext, variant.LocalPhaseSet, accDonType,
                    exonBaseDistance, exonPosition != null ? exonPosition : -1));

            int sampleTotalFrags = 0;
            double samplePsi = PSI_NO_RATE;
            double cohortPsiPerc = PSI_NO_RATE;

            double cohortAvgRate = 0;
            int cohortTotalFrags = 0;

            if(exonPosition != null)
            {
                final int[] spliceSiteData = mSpliceSiteCache.getSampleSpliceSiteData(variant.Chromosome, exonPosition);

                if(spliceSiteData != null)
                {
                    sampleTotalFrags = spliceSiteData[SS_SUPPORT] + spliceSiteData[SS_TRAVERSED];
                    samplePsi = calcSupportRate(spliceSiteData);
                    cohortPsiPerc = mSpliceSiteCache.getCohortSpliceSitePsiPercentile(variant.Chromosome, exonPosition, samplePsi);
                }

                final int[] cohortData = mSpliceSiteCache.getCohortSpliceSiteData(variant.Chromosome, exonPosition);

                if(cohortData != null)
                {
                    cohortAvgRate = calcSupportRate(cohortData);
                    cohortTotalFrags = cohortData[SS_SUPPORT] + cohortData[SS_TRAVERSED];
                }
            }

            if(altSJ != null)
            {
                mWriter.write(String.format(",%s,%s,%s,%d,%d,%d,%d,%d",
                        altSJ.Type, altSJ.RegionContexts[SE_START], altSJ.RegionContexts[SE_END],
                        altSJ.SpliceJunction[SE_START], altSJ.SpliceJunction[SE_END],
                        altSJ.FragmentCount, altSJ.DepthCounts[SE_START], altSJ.DepthCounts[SE_END]));

                mWriter.write(String.format(",%s,%s,%d,%s",
                        getBaseContext(variant, altSJ.SpliceJunction[SE_START], geneData.Strand == POS_STRAND ? DONOR : ACCEPTOR),
                        getBaseContext(variant, altSJ.SpliceJunction[SE_END], geneData.Strand == POS_STRAND ? ACCEPTOR : DONOR),
                        mDataCache.getCohortAltSjFrequency(altSJ), transDataStr));

                mWriter.write(String.format(",%.4f,%d,%.4f,%.4f,%d",
                        samplePsi, sampleTotalFrags, cohortPsiPerc, cohortAvgRate, cohortTotalFrags));
            }
            else
            {
                mWriter.write(String.format(",%s,%s,%s,-1,-1,0,0,0,,,0,,%.4f,%d,%.4f,%.4f,%d",
                        UNKNOWN, AltSpliceJunctionContext.UNKNOWN, AltSpliceJunctionContext.UNKNOWN,
                        samplePsi, sampleTotalFrags, cohortPsiPerc, cohortAvgRate, cohortTotalFrags));
            }

            mWriter.newLine();
        }
        catch(IOException e)
        {
            ISF_LOGGER.error("failed to write splice variant matching file: {}", e.toString());
        }
    }

}
