package com.hartwig.hmftools.esvee.utils;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.deriveRefGenomeVersion;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_ORIENTATION;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.esvee.AssemblyConfig.BWA_LIB_PATH;
import static com.hartwig.hmftools.esvee.AssemblyConfig.SV_LOGGER;
import static com.hartwig.hmftools.esvee.common.FileCommon.APP_NAME;
import static com.hartwig.hmftools.esvee.common.FileCommon.REF_GENOME_IMAGE_EXTENSION;
import static com.hartwig.hmftools.esvee.alignment.AlignmentWriter.FLD_ALIGNED_BASES;
import static com.hartwig.hmftools.esvee.alignment.AlignmentWriter.FLD_CIGAR;
import static com.hartwig.hmftools.esvee.alignment.AlignmentWriter.FLD_FLAGS;
import static com.hartwig.hmftools.esvee.alignment.AlignmentWriter.FLD_MAP_QUAL;
import static com.hartwig.hmftools.esvee.alignment.AlignmentWriter.FLD_MD_TAG;
import static com.hartwig.hmftools.esvee.alignment.AlignmentWriter.FLD_NMATCHES;
import static com.hartwig.hmftools.esvee.alignment.AlignmentWriter.FLD_REF_LOCATION;
import static com.hartwig.hmftools.esvee.alignment.AlignmentWriter.FLD_SCORE;
import static com.hartwig.hmftools.esvee.alignment.AlignmentWriter.FLD_SEQUENCE_COORDS;
import static com.hartwig.hmftools.esvee.alignment.AlignmentWriter.FLD_XA_TAG;
import static com.hartwig.hmftools.esvee.alignment.BwaAligner.loadAlignerLibrary;

import static htsjdk.samtools.CigarOperator.M;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.utils.config.ConfigUtils;
import com.hartwig.hmftools.esvee.alignment.AlignData;
import com.hartwig.hmftools.esvee.alignment.BwaAligner;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;

public class BwaTester
{
    private static final String INPUT_FILE = "input_file";
    private static final String OUTPUT_FILE = "output_file";

    private final String mInputFile;
    private final String mOutputFile;

    private final RefGenomeSource mRefGenome;
    private final RefGenomeVersion mRefGenomeVersion;
    private final BwaAligner mAligner;

    public BwaTester(final ConfigBuilder configBuilder)
    {
        String refGenomeFile = configBuilder.getValue(REF_GENOME);
        mRefGenome = loadRefGenome(refGenomeFile);
        mRefGenomeVersion = deriveRefGenomeVersion(mRefGenome);

        mInputFile = configBuilder.getValue(INPUT_FILE);
        mOutputFile = configBuilder.getValue(OUTPUT_FILE);

        String bwaLibPath = configBuilder.getValue(BWA_LIB_PATH);

        loadAlignerLibrary(bwaLibPath);

        String refGenomeImageFile = refGenomeFile + REF_GENOME_IMAGE_EXTENSION;
        mAligner = new BwaAligner(refGenomeImageFile);
    }

    public void run()
    {
        if(mInputFile == null || mOutputFile == null || !Files.exists(Paths.get(mInputFile)))
        {
            SV_LOGGER.error("invalid input file");
            System.exit(1);
        }

        try
        {
            BufferedWriter writer = createBufferedWriter(mOutputFile);

            StringJoiner sj = new StringJoiner(TSV_DELIM);
            sj.add(FLD_CHROMOSOME).add(FLD_POS_START).add(FLD_POS_END);

            sj.add(FLD_REF_LOCATION);
            sj.add(FLD_SEQUENCE_COORDS);
            sj.add(FLD_MAP_QUAL);
            sj.add(FLD_CIGAR);
            sj.add(FLD_ORIENTATION);
            sj.add(FLD_ALIGNED_BASES);
            sj.add(FLD_SCORE);
            sj.add(FLD_FLAGS);
            sj.add(FLD_NMATCHES);
            sj.add(FLD_XA_TAG);
            sj.add(FLD_MD_TAG);

            writer.write(sj.toString());
            writer.newLine();

            List<String> lines = Files.readAllLines(Paths.get(mInputFile));
            lines.remove(0);

            SV_LOGGER.info("test BWA alignment tests for {} regions", lines.size());

            int chrIndex = 0;
            int posStartIndex = 1;
            int posEndIndex = 2;
            Integer sequenceIndex = 3;

            for(String line : lines)
            {
                String[] values = line.split(TSV_DELIM);

                ChrBaseRegion region = new ChrBaseRegion(values[chrIndex], Integer.parseInt(values[posStartIndex]), Integer.parseInt(values[posEndIndex]));
                String sequence = values.length > sequenceIndex ? values[sequenceIndex] : null;

                testRegion(region, sequence, writer);
            }

            writer.close();
        }
        catch(Exception e)
        {
            SV_LOGGER.error("error reading / writing: {}", e.toString());
            System.exit(1);
        }

        SV_LOGGER.info("BWA alignment tests complete");
    }

    private void testRegion(final ChrBaseRegion region, @Nullable final String sequence, final BufferedWriter writer) throws IOException
    {
        String testSequence = sequence != null ? sequence : mRefGenome.getBaseString(region.Chromosome, region.start(), region.end());

        SV_LOGGER.debug("testing region({})", region);

        List<BwaMemAlignment> bwaAlignments = mAligner.alignSequence(testSequence.getBytes());

        SV_LOGGER.debug("region({}) found {} alignments", region, bwaAlignments.size());

        List<AlignData> alignments = bwaAlignments.stream()
                .map(x -> AlignData.from(x, mRefGenomeVersion))
                .filter(x -> x != null).collect(Collectors.toList());

        for(AlignData alignment : alignments)
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(region.Chromosome);
            sj.add(String.valueOf(region.start()));
            sj.add(String.valueOf(region.end()));

            sj.add(alignment.RefLocation.toString());
            sj.add(format("%d-%d", alignment.rawSequenceStart(), alignment.rawSequenceEnd()));
            sj.add(String.valueOf(alignment.MapQual));
            sj.add(String.valueOf(alignment.Cigar));
            sj.add(String.valueOf(alignment.orientation()));

            Cigar cigar = CigarUtils.cigarFromStr(alignment.Cigar);
            int alignedBases = cigar.getCigarElements().stream().filter(x -> x.getOperator() == M).mapToInt(x -> x.getLength()).sum();
            sj.add(String.valueOf(alignedBases));

            sj.add(String.valueOf(alignment.Score));
            sj.add(String.valueOf(alignment.Flags));
            sj.add(String.valueOf(alignment.NMatches));
            sj.add(alignment.XaTag);
            sj.add(alignment.MdTag);

            writer.write(sj.toString());
            writer.newLine();
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);

        configBuilder.addPath(INPUT_FILE, true, "Input regions and sequences");
        configBuilder.addConfigItem(OUTPUT_FILE, true, "Output file with alignment results");
        configBuilder.addPath(BWA_LIB_PATH, false, "Path to BWA library");

        ConfigUtils.addLoggingOptions(configBuilder);
        addRefGenomeFile(configBuilder, true);

        configBuilder.checkAndParseCommandLine(args);

        BwaTester bwaTester = new BwaTester(configBuilder);

        bwaTester.run();
    }
}
