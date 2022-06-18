package com.hartwig.hmftools.neo.utils;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.addEnsemblDir;
import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_0;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_2;
import static com.hartwig.hmftools.common.gene.CodingBaseData.PHASE_NONE;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME_CFG_DESC;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.ConfigUtils.loadGeneIdsFile;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.neo.NeoCommon.NE_LOGGER;
import static com.hartwig.hmftools.neo.bind.BindConstants.MIN_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.BindConstants.REF_PEPTIDE_LENGTH;
import static com.hartwig.hmftools.neo.bind.FlankCounts.FLANK_AA_COUNT;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Set;

import com.beust.jcommander.internal.Sets;
import com.hartwig.hmftools.common.codon.AminoAcids;
import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptAminoAcids;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

public class MissensePeptideWriter
{
    private final EnsemblDataCache mEnsemblDataCache;
    private final RefGenomeInterface mRefGenome;
    private final List<String> mGeneIds;

    private final int[] mPeptideLengths;
    private final int mFlankLength;
    private final Set<String> mUniqiuePeptides;

    private BufferedWriter mWriter;

    private static final String GENE_ID_FILE = "gene_id_file";
    private static final String OUTPUT_FILE = "output_file";

    public MissensePeptideWriter(final CommandLine cmd)
    {
        mEnsemblDataCache = new EnsemblDataCache(cmd.getOptionValue(ENSEMBL_DATA_DIR), RefGenomeVersion.V37);
        mEnsemblDataCache.setRequiredData(true, false, false, true);

        mRefGenome = loadRefGenome(cmd.getOptionValue(REF_GENOME));

        mPeptideLengths = new int[] { MIN_PEPTIDE_LENGTH, REF_PEPTIDE_LENGTH };
        mFlankLength = FLANK_AA_COUNT;
        mUniqiuePeptides = Sets.newHashSet();

        mWriter = initialiseWriter(cmd.getOptionValue(OUTPUT_FILE));

        mGeneIds = loadGeneIdsFile(cmd.getOptionValue(GENE_ID_FILE));
    }

    public void run()
    {
        if(mGeneIds.isEmpty())
        {
            NE_LOGGER.error("no gene IDs specified");
            System.exit(1);
        }

        mEnsemblDataCache.load(true);
        mEnsemblDataCache.loadTranscriptData(mGeneIds);

        NE_LOGGER.info("writing missense peptides for {} genes", mGeneIds.size());

        int geneCount = 0;

        for(String geneId : mGeneIds)
        {
            GeneData geneData = mEnsemblDataCache.getGeneDataById(geneId);
            List<TranscriptData> transDataList = mEnsemblDataCache.getTranscripts(geneData.GeneId);

            if(transDataList == null)
                continue;

            for(TranscriptData transData : transDataList)
            {
                if(!transData.IsCanonical)
                    continue;

                mUniqiuePeptides.clear();
                processTranscript(geneData, transData);
            }

            ++geneCount;
        }

        closeBufferedWriter(mWriter);

        NE_LOGGER.info("wrote {} gene missense peptide sets", geneCount);
    }

    private void processTranscript(final GeneData geneData, final TranscriptData transData)
    {
        if(transData.CodingStart == null)
            return;

        boolean inCoding = false;

        if(transData.Strand == POS_STRAND)
        {
            StringBuilder codingBases = new StringBuilder();

            for(int i = 0; i < transData.exons().size(); ++i)
            {
                ExonData exon = transData.exons().get(i);

                if(exon.End < transData.CodingStart)
                    continue;

                if(exon.Start > transData.CodingEnd)
                    break;

                int exonCodingStart = max(transData.CodingStart, exon.Start);
                int exonCodingEnd = min(transData.CodingEnd, exon.End);

                if(!inCoding)
                {
                    inCoding = true;

                    if(exon.Start == transData.CodingStart)
                    {
                        int startPhase = exon.PhaseStart == PHASE_NONE ? PHASE_0 : exon.PhaseStart;

                        if(startPhase == PHASE_2)
                            exonCodingStart += 2;
                        else if(startPhase == PHASE_0)
                            exonCodingStart += 1;
                    }
                }

                codingBases.append(mRefGenome.getBaseString(geneData.Chromosome, exonCodingStart, exonCodingEnd));
            }

            produceMissensePeptides(geneData, transData, codingBases.toString());
        }
        else
        {
            String codingBases = "";
            for(int i = transData.exons().size() - 1; i >= 0; --i)
            {
                ExonData exon = transData.exons().get(i);

                if(exon.Start > transData.CodingEnd)
                    continue;

                if(exon.End < transData.CodingStart)
                    break;

                int exonCodingStart = max(transData.CodingStart, exon.Start);
                int exonCodingEnd = min(transData.CodingEnd, exon.End);

                if(!inCoding)
                {
                    inCoding = true;

                    if(exon.End == transData.CodingEnd)
                    {
                        int startPhase = exon.PhaseStart == PHASE_NONE ? PHASE_0 : exon.PhaseStart;

                        if(startPhase == PHASE_2)
                            exonCodingStart -= 2;
                        else if(startPhase == PHASE_0)
                            exonCodingStart -= 1;
                    }
                }

                codingBases = mRefGenome.getBaseString(geneData.Chromosome, exonCodingStart, exonCodingEnd) + codingBases;
            }

            produceMissensePeptides(geneData, transData, Nucleotides.reverseStrandBases(codingBases));
        }
    }

    private void produceMissensePeptides(
            final GeneData geneData, final TranscriptData transData, final String codingBases)
    {
        int codonCount = codingBases.length() / 3;

        // create peptides for the specified length range and with flanks
        for(int codonIndex = mFlankLength; codonIndex < codonCount - mFlankLength; ++codonIndex)
        {
            int codonStartBaseIndex = codonIndex * 3;

            // cycle through each missense variant
            String refCodon = codingBases.substring(codonStartBaseIndex, codonStartBaseIndex + 3);
            String codonRefAminoAcid = AminoAcids.findAminoAcidForCodon(refCodon);

            for(int codonBaseIndex = 0; codonBaseIndex <= 2; ++codonBaseIndex)
            {
                char codonRefBase = refCodon.charAt(codonBaseIndex);

                for(char dnaBase : Nucleotides.DNA_BASES)
                {
                    if(dnaBase == codonRefBase)
                        continue;

                    String altCodon = "";

                    for(int i = 0; i <= 2; ++i)
                    {
                        if(i == codonBaseIndex)
                            altCodon += dnaBase;
                        else
                            altCodon += codingBases.charAt(codonIndex * 3 + i);
                    }

                    String codonAltAminoAcid = AminoAcids.findAminoAcidForCodon(altCodon);

                    if(codonAltAminoAcid == null || codonRefAminoAcid.equals(codonAltAminoAcid))
                        continue;

                    generatePeptides(geneData, transData, codingBases, codonIndex, codonAltAminoAcid, refCodon, altCodon);
                }
            }
        }
    }

    private void generatePeptides(
            final GeneData geneData, final TranscriptData transData, final String codingBases, int codonIndex, final String altAminoAcid,
            final String refCodon, final String altCodon)
    {
        for(int pepLen = mPeptideLengths[0]; pepLen <= mPeptideLengths[1]; ++pepLen)
        {
            // if say length is 8 and the current codon index is 15 with flank of 3, then the peptide + flanks needs to cover
            // the codon range 8 -> 15 to 15 -> 22, which flanks becomes 5 -> 15 up to 15 -> 25
            int peptideCodonStartIndex = max(0, codonIndex - pepLen - mFlankLength + 1);

            NE_LOGGER.trace("gene({}) pepLen({}) codonIndex({}) codonRangeStart({})",
                    geneData.GeneName, pepLen, codonIndex, peptideCodonStartIndex);

            for(int pepStartCodonIndex = peptideCodonStartIndex; pepStartCodonIndex <= codonIndex - mFlankLength; ++pepStartCodonIndex)
            {
                int pepPlusFlanksLen = pepLen + mFlankLength * 2;

                int peptideEndCodingBaseIndex = (pepStartCodonIndex + pepPlusFlanksLen - 1) * 3;
                if(peptideEndCodingBaseIndex >= codingBases.length())
                    break;

                String peptide = "";

                for(int pepIndex = 0; pepIndex < pepPlusFlanksLen; ++pepIndex)
                {
                    int codingBaseCodonIndex = pepStartCodonIndex + pepIndex;

                    if(codingBaseCodonIndex == codonIndex)
                    {
                        peptide += altAminoAcid;
                    }
                    else
                    {
                        int pepCodonBaseIndex = (pepStartCodonIndex + pepIndex) * 3;
                        String codon = codingBases.substring(pepCodonBaseIndex, pepCodonBaseIndex + 3);
                        String aminoAcid = AminoAcids.findAminoAcidForCodon(codon);;

                        if(aminoAcid == null)
                        {
                            peptide = "";
                            break;
                        }

                        peptide += aminoAcid;
                    }
                }

                if(peptide.isEmpty())
                    continue;

                NE_LOGGER.trace("gene({}) pepLen({}) codonIndex({}) codonRange({} -> {}) peptide({})",
                        geneData.GeneName, pepLen, codonIndex, pepStartCodonIndex, pepStartCodonIndex + pepPlusFlanksLen - 1, peptide);

                writePeptideData(geneData.GeneId, geneData.GeneName, transData.TransName, codonIndex, peptide, refCodon, altCodon);
            }
        }
    }

    private static BufferedWriter initialiseWriter(final String outputFile)
    {
        NE_LOGGER.info("writing results to {}", outputFile);

        try
        {
            BufferedWriter writer = createBufferedWriter(outputFile, false);

            writer.write("GeneId,GeneName,TransName,Peptide,UpFlank,DownFlank,CodonIndex,RefCodon,AltCodon");
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to initialise CSV file output: {}", e.toString());
            return null;
        }
    }

    private void writePeptideData(
            final String geneId, final String geneName, final String transName, int codonIndex, final String peptidePlusFlanks,
            final String refCodon, final String altCodon)
    {
        if(peptidePlusFlanks.length() < mFlankLength * 2 + 1)
            return;

        String peptide = peptidePlusFlanks.substring(mFlankLength, peptidePlusFlanks.length() - mFlankLength);

        if(mUniqiuePeptides.contains(peptide))
            return;

        mUniqiuePeptides.add(peptide);

        String upFlank = peptidePlusFlanks.substring(0, mFlankLength);
        String downFlank = peptidePlusFlanks.substring(peptidePlusFlanks.length() - mFlankLength);

        try
        {
            mWriter.write(String.format("%s,%s,%s,%s,%s,%s,%d,%s,%s",
                    geneId, geneName, transName, peptide, upFlank, downFlank, codonIndex, refCodon, altCodon));
            mWriter.newLine();
        }
        catch(IOException e)
        {
            NE_LOGGER.error("failed to write peptide data: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = new Options();
        addEnsemblDir(options);
        options.addOption(REF_GENOME, true, REF_GENOME_CFG_DESC);
        options.addOption(GENE_ID_FILE, true, "Gene IDs file");
        options.addOption(OUTPUT_FILE, true, "Output filename");
        options.addOption(OUTPUT_DIR, true, "Output directory");
        addLoggingOptions(options);

        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        MissensePeptideWriter missensePeptideWriter = new MissensePeptideWriter(cmd);
        missensePeptideWriter.run();
    }

    @NotNull
    public static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
