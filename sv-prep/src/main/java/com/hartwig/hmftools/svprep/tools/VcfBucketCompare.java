package com.hartwig.hmftools.svprep.tools;

import static java.lang.Math.abs;
import static java.lang.String.format;

import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.IHOMPOS;
import static com.hartwig.hmftools.common.sv.StructuralVariantFactory.createSingleBreakend;
import static com.hartwig.hmftools.common.sv.StructuralVariantType.SGL;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_ID;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.addOutputOptions;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.parseOutputDir;
import static com.hartwig.hmftools.common.utils.sv.BaseRegion.positionWithin;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_CHROMOSOMES;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_CHROMOSOMES_DESC;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificChromsomes;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.svprep.SvCommon.DELIM;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConfig.SAMPLE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.vcf.VCFCodec;

public class VcfBucketCompare
{
    private final String mSampleId;
    private final String mVcfFilename;
    private final String mJunctionsFilename;
    private final String mOutputDir;
    private final String mOutputId;
    private final List<String> mSpecificChromosomes;

    private final Map<String,List<JunctionData>> mChrJunctions;
    private final Map<String,List<RemoteJunction>> mChrRemoteJunctions;
    private StructuralVariantFactory mSvFactory;
    private final Map<MatchType,Integer> mMatchCounts;

    private int mVcfBreakends;
    private final BufferedWriter mWriter;

    private static final String VCF_FILE = "vcf_file";
    private static final String BUCKET_FILE = "bucket_file";
    private static final String JUNCTION_FILE = "junction_file";

    public VcfBucketCompare(final CommandLine cmd)
    {
        mSampleId = cmd.getOptionValue(SAMPLE);
        mVcfFilename = cmd.getOptionValue(VCF_FILE);
        mJunctionsFilename = cmd.getOptionValue(JUNCTION_FILE);
        mOutputDir = parseOutputDir(cmd);
        mOutputId = cmd.getOptionValue(OUTPUT_ID);

        mChrJunctions = Maps.newHashMap();
        mChrRemoteJunctions = Maps.newHashMap();
        mMatchCounts = Maps.newHashMap();
        mSvFactory = new StructuralVariantFactory(new CompoundFilter(false));
        mWriter = initialiseWriter();
        mVcfBreakends = 0;

        mSpecificChromosomes = loadSpecificChromsomes(cmd);
    }

    public void run()
    {
        if(mVcfFilename == null)
        {
            SV_LOGGER.error("missing VCF");
            System.exit(1);
        }

        if(!loadJunctionData(mJunctionsFilename))
        {
            SV_LOGGER.error("invalid junctions data file");
            System.exit(1);
        }

        SV_LOGGER.info("comparing junctions data to VCF({})", mVcfFilename);

        final AbstractFeatureReader<VariantContext, LineIterator> reader = AbstractFeatureReader.getFeatureReader(
                mVcfFilename, new VCFCodec(), false);

        // VCFHeader vcfHeader = (VCFHeader)reader.getHeader();

        try
        {
            for(VariantContext variantContext : reader.iterator())
            {
                boolean isSgl = StructuralVariantFactory.isSingleBreakend(variantContext);

                StructuralVariant sv = null;

                if(isSgl)
                {
                    sv = createSingleBreakend(variantContext);
                }
                else
                {
                    int currentSvCount = mSvFactory.results().size();
                    mSvFactory.addVariantContext(variantContext);

                    // check if both breakends have now been encountered
                    if(currentSvCount == mSvFactory.results().size())
                        continue;

                    if(mSvFactory.results().isEmpty())
                        continue;

                    sv = mSvFactory.results().get(0);
                    mSvFactory.results().remove(0);
                }

                processBreakend(sv, variantContext);
            }
        }
        catch(IOException e)
        {
            SV_LOGGER.error("error reading vcf({}): {}", mVcfFilename, e.toString());
            System.exit(1);
        }

        StringJoiner matches = new StringJoiner(", ");
        mMatchCounts.entrySet().forEach(x -> matches.add(format("%s=%d", x.getKey(), x.getValue())));
        SV_LOGGER.info("compared {} breakends, matchCounts: {}", mVcfBreakends, matches.toString());

        closeBufferedWriter(mWriter);
    }

    private enum MatchType
    {
        NONE,
        SV,
        JUNCTION,
        NEAR_JUNCTION,
        REMOTE;
    }

    private void addMatchType(MatchType matchType)
    {
        Integer count = mMatchCounts.get(matchType);
        mMatchCounts.put(matchType, count != null ? count + 1 : 1);
    }

    private static int MAX_JUNCTION_DISTANCE = 1000;
    private static int REMOTE_MATCH_BUFFER = 5;

    private void processBreakend(final StructuralVariant sv, final VariantContext variantContext)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(se == SE_END && sv.type() == SGL)
                continue;

            StructuralVariantLeg svLeg = se == SE_START ? sv.start() : sv.end();
            StructuralVariantLeg svOtherLeg = se == SE_START ? sv.end() : sv.start();

            if(!mSpecificChromosomes.isEmpty() && !mSpecificChromosomes.contains(svLeg.chromosome()))
                continue;

            ++mVcfBreakends;

            int[] homology = {0, 0};

            if(variantContext.hasAttribute(IHOMPOS))
            {
                final List<Integer> ihompos = variantContext.getAttributeAsIntList(IHOMPOS, 0);
                homology[0] = ihompos.get(0);
                homology[1] = ihompos.get(1);
            }

            int svPosStart = svLeg.position() + homology[0];
            int svPosEnd = svLeg.position() + homology[1];

            List<JunctionData> junctions = mChrJunctions.get(svLeg.chromosome());

            MatchType matchType = MatchType.NONE;
            JunctionData nearestJunction = null;
            RemoteJunction matchedRemoteJunction = null;
            int nearestDistance = -1;

            for(JunctionData junctionData : junctions)
            {
                if(junctionData.matches(svPosStart, svPosEnd, svLeg.orientation()))
                {
                    matchType = MatchType.JUNCTION;
                    nearestJunction = junctionData;

                    // check for a remote match

                    if(svOtherLeg != null)
                    {
                        matchedRemoteJunction = junctionData.findRemoteMatch(
                                svOtherLeg.chromosome(), svOtherLeg.position(), svOtherLeg.orientation());
                        matchType = MatchType.SV;
                    }

                    break;
                }

                if(junctionData.Orientation != svLeg.orientation())
                    continue;

                int posDiff = abs(svLeg.position() - junctionData.Position);
                if(posDiff > MAX_JUNCTION_DISTANCE)
                    continue;

                if(nearestJunction == null || posDiff < nearestDistance)
                {
                    nearestJunction = junctionData;
                    nearestDistance = posDiff;
                    matchType = MatchType.NEAR_JUNCTION;
                }
            }

            if((nearestJunction == null || matchType == MatchType.NEAR_JUNCTION) && svOtherLeg != null)
            {
                // check the remotes
                List<RemoteJunction> remoteJunctions = mChrRemoteJunctions.get(svOtherLeg.chromosome());

                if(remoteJunctions != null)
                {
                    RemoteJunction remoteJunction = remoteJunctions.stream()
                            .filter(x -> x.matches(svOtherLeg.chromosome(), svOtherLeg.position(), svOtherLeg.orientation()))
                            .findFirst().orElse(null);

                    if(remoteJunction != null && remoteJunction.Junction.Chromosome.equals(svLeg.chromosome())
                    && remoteJunction.Junction.matches(svPosStart, svPosEnd, svLeg.orientation()))
                    {
                        matchedRemoteJunction = remoteJunction;
                        matchType = MatchType.REMOTE;
                        nearestJunction = remoteJunction.Junction;
                    }
                }
            }

            writeMatchData(sv, svLeg, matchType, nearestJunction, matchedRemoteJunction);

            addMatchType(matchType);
        }
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String fileName = mOutputDir + mSampleId + ".sv_prep";

            if(mOutputId != null)
                fileName += "." + mOutputId;

            fileName += ".junction_vcf_compare.csv";

            SV_LOGGER.info("writing comparison file: {}", fileName);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("VcfId,Type,Chromosome,Position,Orientation,InsSeqLength,MateChromosome,MatePosition");
            writer.write(",Filter,Qual,NormalFrags,TumorFrags");
            writer.write(",MatchType,JunctionPos,JunctionFrags,JunctionSupportReads");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private void writeMatchData(
            final StructuralVariant sv, final StructuralVariantLeg svLeg, MatchType matchType,
            final JunctionData junctionData, final RemoteJunction remoteJunction)
    {
        try
        {
            final StructuralVariantLeg otherLeg = sv.type() != null ? (sv.start() == svLeg ? sv.end() : sv.start()) : null;

            mWriter.write(format("%s,%s,%s,%d,%d,%d,%s,%d",
                    sv.id(), sv.type(), svLeg.chromosome(), svLeg.position(), svLeg.orientation(), sv.insertSequence().length(),
                    otherLeg != null ? otherLeg.chromosome() : "", otherLeg != null ? otherLeg.position() : -1));

            mWriter.write(format(",%s,%.0f,%d,%d",
                    sv.filter(), sv.qualityScore(), svLeg.normalVariantFragmentCount(), svLeg.tumorVariantFragmentCount()));

            if(junctionData == null)
            {
                mWriter.write(format(",%s,0,0,0", MatchType.NONE));
            }
            else
            {
                mWriter.write(format(",%s,%d,%d,%d", matchType, junctionData.Position, junctionData.ExactFragments, junctionData.SupportReads));
            }

            mWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write data: {}", e.toString());
        }
    }

    private boolean loadJunctionData(final String filename)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(line, DELIM);

            // Chromosome,BucketStart,BucketEnd,Position,Orientation,Fragments,SupportingReads,Hotspot,InitialReadId,RemoteJunctionCount,RemoteChromosome,RemotePosition,RemoteOrientation
            // 1,118516001,118517000,118516222,-1,30,6,false,A00624:8:HHKYHDSXX:4:2243:4562:28494,2,1,118516186,1

            int chrIndex = fieldsIndexMap.get("Chromosome");
            int posIndex = fieldsIndexMap.get("Position");
            int orientIndex = fieldsIndexMap.get("Orientation");
            int fragsIndex = fieldsIndexMap.get("Fragments");
            int supportIndex = fieldsIndexMap.get("SupportingReads");
            int remoteChrIndex = fieldsIndexMap.get("RemoteChromosome");
            int remotePosIndex = fieldsIndexMap.get("RemotePosition");
            int remoteOrientIndex = fieldsIndexMap.get("RemoteOrientation");

            int junctionCount = 0;
            JunctionData junctionData = null;

            while ((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(DELIM, -1);

                String chromosome = values[chrIndex];

                List<JunctionData> junctions = mChrJunctions.get(chromosome);

                if(junctions == null)
                {
                    junctions = Lists.newArrayList();
                    mChrJunctions.put(chromosome, junctions);
                }

                int position = Integer.parseInt(values[posIndex]);
                byte orientation = Byte.parseByte(values[orientIndex]);

                if(junctionData == null || junctionData.Position != position || junctionData.Orientation != orientation)
                {
                    junctionData = new JunctionData(
                            chromosome, position, orientation, Integer.parseInt(values[fragsIndex]), Integer.parseInt(values[supportIndex]));
                    junctions.add(junctionData);
                    ++junctionCount;
                }

                String remoteChr = values[remoteChrIndex];

                if(!remoteChr.isEmpty())
                {
                    RemoteJunction remoteJunction = new RemoteJunction(
                            junctionData, remoteChr, Integer.parseInt(values[remotePosIndex]), Byte.parseByte(values[remoteOrientIndex]));

                    junctionData.RemoteJunctions.add(remoteJunction);

                    List<RemoteJunction> remoteJunctions = mChrRemoteJunctions.get(remoteChr);

                    if(remoteJunctions == null)
                    {
                        remoteJunctions = Lists.newArrayList();
                        mChrRemoteJunctions.put(remoteChr, remoteJunctions);
                    }

                    if(remoteJunctions.stream().noneMatch(x -> x.matches(remoteJunction)))
                    {
                        remoteJunctions.add(remoteJunction);
                    }
                }
            }

            SV_LOGGER.debug("loaded {} junction data records from file: {}", junctionCount, filename);
        }
        catch(IOException exception)
        {
            SV_LOGGER.error("failed to read junction data file({})", filename, exception.toString());
            return false;
        }

        return true;
    }

    // Chromosome,BucketStart,BucketEnd,Position,Orientation,Fragments,SupportingReads,Hotspot,InitialReadId,RemoteJunctionCount,RemoteChromosome,RemotePosition,RemoteOrientation
    //  1,118516001,118517000,118516222,-1,30,6,false,A00624:8:HHKYHDSXX:4:2243:4562:28494,2,1,118516186,1

    private class JunctionData
    {
        public final String Chromosome;
        public final int Position;
        public final byte Orientation;
        public final int ExactFragments;
        public final int SupportReads;
        public final List<RemoteJunction> RemoteJunctions;

        public JunctionData(final String chromosome, final int position, final byte orientation, final int exactFragments, final int supportReads)
        {
            Chromosome = chromosome;
            Position = position;
            Orientation = orientation;
            ExactFragments = exactFragments;
            SupportReads = supportReads;
            RemoteJunctions = Lists.newArrayList();
        }

        public boolean matches(int posStart, int posEnd, byte orientation)
        {
            return Orientation == orientation && positionWithin(Position, posStart, posEnd);
        }

        public RemoteJunction findRemoteMatch(final String chromosome, final int position, final byte orientation)
        {
            return RemoteJunctions.stream().filter(x -> x.matches(chromosome, position, orientation)).findFirst().orElse(null);
        }

        public String toString()
        {
            return format("%s:%d %d) frags(exact=%d supp=%d) remotes(%d)",
                    Chromosome, Position, Orientation, ExactFragments, SupportReads, RemoteJunctions.size());
        }

    }

    public class RemoteJunction
    {
        public final JunctionData Junction;
        public final String Chromosome;
        public final int Position;
        public final byte Orientation;
        public int Support;

        public RemoteJunction(final JunctionData junction, final String chromosome, final int position, final byte orientation)
        {
            Chromosome = chromosome;
            Junction = junction;
            Position = position;
            Orientation = orientation;
            Support = 0;
        }

        public boolean matches(final RemoteJunction other)
        {
            return Position == other.Position && Orientation == other.Orientation;
        }

        public boolean matches(final String chromosome, final int position, final byte orientation)
        {
            return Chromosome.equals(chromosome) && abs(Position - position) <= REMOTE_MATCH_BUFFER && Orientation == orientation;
        }

        public String toString() { return format("loc(%s:%d:%d) reads(%d) junction(%s)",
                Chromosome, Position, Orientation, Support, Junction); }
    }


    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Name of the sample");
        options.addOption(VCF_FILE, true, "VCF File");
        options.addOption(BUCKET_FILE, true, "SV prep bucket file");
        options.addOption(JUNCTION_FILE, true, "SV prep bucket file");
        options.addOption(SPECIFIC_CHROMOSOMES, true, SPECIFIC_CHROMOSOMES_DESC);

        addOutputOptions(options);
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        VcfBucketCompare vcfBucketCompare = new VcfBucketCompare(cmd);
        vcfBucketCompare.run();

        SV_LOGGER.info("SV prep buckets vs VCF complete");
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }
}
