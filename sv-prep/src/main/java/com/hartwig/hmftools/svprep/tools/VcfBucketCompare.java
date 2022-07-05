package com.hartwig.hmftools.svprep.tools;

import static java.lang.Math.abs;
import static java.lang.String.format;

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
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_CHROMOSOMES;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.SPECIFIC_CHROMOSOMES_DESC;
import static com.hartwig.hmftools.common.utils.sv.ChrBaseRegion.loadSpecificChromsomes;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_END;
import static com.hartwig.hmftools.common.utils.sv.StartEndIterator.SE_START;
import static com.hartwig.hmftools.svprep.SvCommon.DELIM;
import static com.hartwig.hmftools.svprep.SvCommon.ITEM_DELIM;
import static com.hartwig.hmftools.svprep.SvCommon.SUB_ITEM_DELIM;
import static com.hartwig.hmftools.svprep.SvCommon.SV_LOGGER;
import static com.hartwig.hmftools.svprep.SvConfig.SAMPLE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFactory;
import com.hartwig.hmftools.common.sv.StructuralVariantLeg;
import com.hartwig.hmftools.common.utils.sv.ChrBaseRegion;

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
    private final String mBucketsFilename;
    private final String mOutputDir;
    private final String mOutputId;
    private final List<String> mSpecificChromosomes;

    private final Map<String,List<BucketData>> mChrBuckets;
    private StructuralVariantFactory mSvFactory;

    private int mVcfBreakends;
    private int mMatchedBreakends;
    private final BufferedWriter mWriter;

    private static final String VCF_FILE = "vcf_file";
    private static final String BUCKET_FILE = "bucket_file";

    public VcfBucketCompare(final CommandLine cmd)
    {
        mSampleId = cmd.getOptionValue(SAMPLE);
        mVcfFilename = cmd.getOptionValue(VCF_FILE);
        mBucketsFilename = cmd.getOptionValue(BUCKET_FILE);
        mOutputDir = parseOutputDir(cmd);
        mOutputId = cmd.getOptionValue(OUTPUT_ID);

        mChrBuckets = Maps.newHashMap();
        mSvFactory = new StructuralVariantFactory(new CompoundFilter(false));
        mWriter = initialiseWriter();

        mMatchedBreakends = 0;
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

        if(!loadBucketData(mBucketsFilename))
        {
            SV_LOGGER.error("invalid bucket data file");
            System.exit(1);
        }

        SV_LOGGER.info("comparing bucket data to VCF({})", mVcfFilename);

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

        SV_LOGGER.info("compared {} breakends, matched({})", mVcfBreakends, mMatchedBreakends);

        closeBufferedWriter(mWriter);
    }

    private void processBreakend(final StructuralVariant sv, final VariantContext variantContext)
    {
        for(int se = SE_START; se <= SE_END; ++se)
        {
            if(se == SE_END && sv.type() == SGL)
                continue;

            StructuralVariantLeg svLeg = se == SE_START ? sv.start() : sv.end();

            if(!mSpecificChromosomes.isEmpty() && !mSpecificChromosomes.contains(svLeg.chromosome()))
                continue;

            ++mVcfBreakends;

            List<BucketData> buckets = mChrBuckets.get(svLeg.chromosome());

            if(buckets == null)
            {
                writeMatchData(sv, svLeg, null, null);
                continue;
            }

            BucketData bucketData = buckets.stream().filter(x -> x.Region.containsPosition(svLeg.position())).findFirst().orElse(null);

            if(bucketData == null)
            {
                writeMatchData(sv, svLeg, null, null);
                continue;
            }

            JunctionData junctionData = bucketData.findBreakendMatch(svLeg.position(), svLeg.orientation());
            writeMatchData(sv, svLeg, bucketData, junctionData);
            ++mMatchedBreakends;
        }
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String fileName = mOutputDir + mSampleId + ".sv_prep";

            if(mOutputId != null)
                fileName += "." + mOutputId;

            fileName += ".bucket_vcf_compare.csv";

            SV_LOGGER.info("writing comparison file: {}", fileName);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("VcfId,Type,Chromosome,Position,Orientation,MateChromosome,MatePosition,Qual,NormalFrags,TumorFrags");
            writer.write(",MatchType,BucketJuncFrags,BucketSupportReads");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private enum MatchType
    {
        NONE,
        BUCKET,
        JUNCTION;
    }

    private void writeMatchData(
            final StructuralVariant sv, final StructuralVariantLeg svLeg, final BucketData bucketData, final JunctionData junctionData)
    {
        try
        {
            final StructuralVariantLeg otherLeg = sv.type() != null ? (sv.start() == svLeg ? sv.end() : sv.start()) : null;

            mWriter.write(format("%s,%s,%s,%d,%d,%s,%d",
                    sv.id(), sv.type(), svLeg.chromosome(), svLeg.position(), svLeg.orientation(),
                    otherLeg != null ? otherLeg.chromosome() : "", otherLeg != null ? otherLeg.position() : -1));

            mWriter.write(format(",%.0f,%d,%d",
                    sv.qualityScore(), svLeg.normalVariantFragmentCount(), svLeg.tumorVariantFragmentCount()));

            if(bucketData == null)
            {
                mWriter.write(format(",%s,0,0", MatchType.NONE));
            }
            else if(junctionData == null)
            {
                mWriter.write(format(",%s,0,0", MatchType.BUCKET));
            }
            else
            {
                mWriter.write(format(",%s,%d,%d", MatchType.JUNCTION, junctionData.ExactFragments, junctionData.SupportReads));
            }

            mWriter.newLine();
        }
        catch(IOException e)
        {
            SV_LOGGER.error("failed to write data: {}", e.toString());
        }
    }

    private boolean loadBucketData(final String filename)
    {
        try
        {
            BufferedReader fileReader = new BufferedReader(new FileReader(filename));

            String line = fileReader.readLine();
            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(line, DELIM);

            int chrIndex = fieldsIndexMap.get("Chromosome");
            int posStartIndex = fieldsIndexMap.get("PosStart");
            int posEndIndex = fieldsIndexMap.get("PosEnd");
            int juncDataIndex = fieldsIndexMap.get("Junctions");

            int bucketCount = 0;

            while ((line = fileReader.readLine()) != null)
            {
                final String[] values = line.split(DELIM, -1);

                String chromosome = values[chrIndex];

                List<BucketData> buckets = mChrBuckets.get(chromosome);

                if(buckets == null)
                {
                    buckets = Lists.newArrayList();
                    mChrBuckets.put(chromosome, buckets);
                }

                BucketData bucket = new BucketData(new ChrBaseRegion(
                        chromosome, Integer.parseInt(values[posStartIndex]), Integer.parseInt(values[posEndIndex])));

                String[] junctions = values[juncDataIndex].split(ITEM_DELIM, -1);

                for(String junctionStr : junctions)
                {
                    String[] items = junctionStr.split(SUB_ITEM_DELIM, 4);
                    JunctionData junctionData = new JunctionData(
                            Integer.parseInt(items[0]), Byte.parseByte(items[1]),
                            Integer.parseInt(items[2]), Integer.parseInt(items[3]));
                    bucket.Junctions.add(junctionData);
                }

                buckets.add(bucket);
                ++bucketCount;
            }

            SV_LOGGER.debug("loaded {} bucket data records from file: {}", bucketCount, filename);
        }
        catch(IOException exception)
        {
            SV_LOGGER.error("failed to read bucket data file({})", filename, exception.toString());
            return false;
        }

        return true;
    }

    private static final int MATCH_BUFFER = 5;

    private class BucketData
    {
        public final ChrBaseRegion Region;

        public final List<JunctionData> Junctions;

        public BucketData(final ChrBaseRegion region)
        {
            Region = region;
            Junctions = Lists.newArrayList();
        }

        public JunctionData findBreakendMatch(int position, byte orientation)
        {
            int minDiff = -1;
            JunctionData topMatch = null;

            for(JunctionData junctionData : Junctions)
            {
                if(junctionData.Orientation != orientation)
                    continue;

                int posDiff = abs(junctionData.Position - position);
                if(posDiff > MATCH_BUFFER)
                    continue;

                if(topMatch == null || posDiff < minDiff)
                {
                    topMatch = junctionData;
                    minDiff = posDiff;
                }
            }

            return topMatch;
        }
    }

    private class JunctionData
    {
        public final int Position;
        public final byte Orientation;
        public final int ExactFragments;
        public final int SupportReads;

        public JunctionData(final int position, final byte orientation, final int exactFragments, final int supportReads)
        {
            Position = position;
            Orientation = orientation;
            ExactFragments = exactFragments;
            SupportReads = supportReads;
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        options.addOption(SAMPLE, true, "Name of the sample");
        options.addOption(VCF_FILE, true, "VCF File");
        options.addOption(BUCKET_FILE, true, "SV prep bucket file");
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
