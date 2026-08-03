package com.hartwig.hmftools.panelbuilder.samplevariants.testdata;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.StringJoiner;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.linx.DriverEventType;
import com.hartwig.hmftools.common.linx.ImmutableLinxCluster;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxCluster;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.linx.LinxGermlineDisruption;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.linx.LinxTestFactory;
import com.hartwig.hmftools.common.purple.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.purple.SegmentSupport;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.common.variant.impact.VariantImpactSerialiser;

import htsjdk.samtools.util.BlockCompressedOutputStream;

// Generates fake sample variant input files (Purple + Linx) for PanelBuilder's sample variant feature.
// The data is entirely synthetic: invented variants placed at genome coordinates, with reference bases read from the supplied
// reference genome so that probe generation against that same genome is valid. No patient data is used.
// Loading and filtering do not need the genome; it is only used here to fill reference bases and, when the app runs, to build probe
// sequences. Point the app at the same reference genome that was passed to this generator.
public class SampleVariantsTestData
{
    // A somatic or germline SNV/INDEL. Reference bases come from the genome; the alt is derived per kind.
    public record SmallVariant(
            String chromosome,
            int position,
            Kind kind,
            // For INSERT: inserted bases. For DELETE: number of reference bases deleted (as a count via param). For SNV: unused.
            String insertedBases,
            int deletedBases,
            boolean reported,
            String gene,
            CodingEffect codingEffect,
            GermlineStatus germlineStatus,
            double subclonalLikelihood,
            int repeatCount,
            int tumorRefFragments,
            int tumorAltFragments)
    {
        public enum Kind
        {
            SNV,
            INSERT,
            DELETE
        }

        public static SmallVariant snv(String chromosome, int position, boolean reported, String gene, CodingEffect codingEffect,
                GermlineStatus germlineStatus, double subclonalLikelihood, int refFrags, int altFrags)
        {
            return new SmallVariant(
                    chromosome, position, Kind.SNV, "", 0, reported, gene, codingEffect, germlineStatus,
                    subclonalLikelihood, 0, refFrags, altFrags);
        }

        public static SmallVariant deletion(String chromosome, int position, int deletedBases, boolean reported, String gene,
                CodingEffect codingEffect, GermlineStatus germlineStatus, int refFrags, int altFrags)
        {
            return new SmallVariant(
                    chromosome, position, Kind.DELETE, "", deletedBases, reported, gene, codingEffect, germlineStatus,
                    0, 0, refFrags, altFrags);
        }
    }

    // A somatic structural variant plus the Linx annotation needed to classify it as a specific driver type.
    public record Sv(
            String startChromosome,
            int startPosition,
            byte startOrientation,
            String endChromosome,
            int endPosition,
            byte endOrientation,
            String insertSequence,
            DriverKind driverKind,
            String gene,
            double purpleAf,
            double copyNumberChange,
            double junctionCopyNumber,
            int tumorFragments)
    {
        public enum DriverKind
        {
            FUSION,
            AMP,
            DEL,
            DISRUPTION
        }
    }

    // A germline structural variant (Linx germline disruption + reported breakend).
    public record GermlineSv(
            String startChromosome,
            int startPosition,
            byte startOrientation,
            String endChromosome,
            int endPosition,
            byte endOrientation,
            String insertSequence,
            String gene)
    {
    }

    public record Specs(
            List<SmallVariant> somaticSmallVariants,
            List<SmallVariant> germlineSmallVariants,
            List<Sv> somaticSvs,
            List<GermlineSv> germlineSvs)
    {
    }

    public record OutputDirs(
            String purpleDir,
            String linxDir,
            String linxGermlineDir)
    {
    }

    private static final String REFERENCE_SAMPLE_SUFFIX = "R";

    public static OutputDirs generate(RefGenomeInterface genome, String sampleId, String outputRoot, Specs specs)
    {
        String purpleDir = outputRoot + File.separator + "purple";
        String linxDir = outputRoot + File.separator + "linx";
        String linxGermlineDir = outputRoot + File.separator + "linx_germline";
        createDir(purpleDir);
        createDir(linxDir);
        createDir(linxGermlineDir);

        String refSample = sampleId + REFERENCE_SAMPLE_SUFFIX;

        writeSmallVariantVcf(
                genome, PurpleCommon.purpleSomaticVcfFile(purpleDir, sampleId), refSample, sampleId,
                specs.somaticSmallVariants());
        writeSmallVariantVcf(
                genome, PurpleCommon.purpleGermlineVcfFile(purpleDir, sampleId), refSample, sampleId,
                specs.germlineSmallVariants());
        writeSomaticSvs(genome, purpleDir, linxDir, refSample, sampleId, specs.somaticSvs());
        writeGermlineSvs(genome, linxGermlineDir, sampleId, specs.germlineSvs());

        return new OutputDirs(purpleDir, linxDir, linxGermlineDir);
    }

    private static void writeSmallVariantVcf(RefGenomeInterface genome, String path, String refSample, String tumorSample,
            List<SmallVariant> variants)
    {
        StringBuilder vcf = new StringBuilder();
        vcf.append("##fileformat=VCFv4.2\n");
        vcf.append("##INFO=<ID=REPORTED,Number=0,Type=Flag,Description=\"\">\n");
        vcf.append("##INFO=<ID=PURPLE_GERMLINE,Number=1,Type=String,Description=\"\">\n");
        vcf.append("##INFO=<ID=SUBCL,Number=1,Type=Float,Description=\"\">\n");
        vcf.append("##INFO=<ID=RC_REPC,Number=1,Type=Integer,Description=\"\">\n");
        vcf.append("##INFO=<ID=PURPLE_AF,Number=1,Type=Float,Description=\"\">\n");
        vcf.append("##INFO=<ID=PURPLE_CN,Number=1,Type=Float,Description=\"\">\n");
        vcf.append("##INFO=<ID=PURPLE_VCN,Number=1,Type=Float,Description=\"\">\n");
        vcf.append("##INFO=<ID=IMPACT,Number=.,Type=String,Description=\"\">\n");
        vcf.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">\n");
        vcf.append("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"\">\n");
        vcf.append("##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"\">\n");
        vcf.append("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"\">\n");
        appendContigs(genome, vcf, variants.stream().map(SmallVariant::chromosome).distinct().toList());
        vcf.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t").append(refSample).append('\t').append(tumorSample)
                .append('\n');

        for(SmallVariant variant : variants)
        {
            String[] refAlt = resolveRefAlt(genome, variant);
            String ref = refAlt[0];
            String alt = refAlt[1];

            int depth = variant.tumorRefFragments() + variant.tumorAltFragments();
            double af = depth == 0 ? 0 : (double) variant.tumorAltFragments() / depth;

            StringJoiner info = new StringJoiner(";");
            if(variant.reported())
            {
                info.add("REPORTED");
            }
            info.add("PURPLE_GERMLINE=" + variant.germlineStatus().toString());
            info.add(format("SUBCL=%.3f", variant.subclonalLikelihood()));
            info.add("RC_REPC=" + variant.repeatCount());
            info.add(format("PURPLE_AF=%.3f", af));
            info.add("PURPLE_CN=2.0");
            info.add("PURPLE_VCN=1.0");
            info.add("IMPACT=" + impactString(variant.gene(), variant.codingEffect()));

            String refGenotype = "0/0:20,0:0.0:20";
            String tumorGenotype = format(
                    "0/1:%d,%d:%.3f:%d",
                    variant.tumorRefFragments(), variant.tumorAltFragments(), af, depth);

            vcf.append(variant.chromosome()).append('\t')
                    .append(variant.position()).append('\t')
                    .append('.').append('\t')
                    .append(ref).append('\t')
                    .append(alt).append('\t')
                    .append("100").append('\t')
                    .append("PASS").append('\t')
                    .append(info).append('\t')
                    .append("GT:AD:AF:DP").append('\t')
                    .append(refGenotype).append('\t')
                    .append(tumorGenotype).append('\n');
        }

        writeBgzip(path, vcf.toString());
    }

    // Resolves the reference and alt allele strings, reading actual reference bases from the genome so the variant is valid against it.
    private static String[] resolveRefAlt(RefGenomeInterface genome, SmallVariant variant)
    {
        switch(variant.kind())
        {
            case SNV:
            {
                String ref = genome.getBaseString(variant.chromosome(), variant.position(), variant.position()).toUpperCase();
                String alt = differentBase(ref);
                return new String[] { ref, alt };
            }
            case INSERT:
            {
                String anchor = genome.getBaseString(variant.chromosome(), variant.position(), variant.position()).toUpperCase();
                return new String[] { anchor, anchor + variant.insertedBases() };
            }
            case DELETE:
            {
                int endPos = variant.position() + variant.deletedBases();
                String ref = genome.getBaseString(variant.chromosome(), variant.position(), endPos).toUpperCase();
                String alt = ref.substring(0, 1);
                return new String[] { ref, alt };
            }
            default:
                throw new IllegalArgumentException();
        }
    }

    private static void writeSomaticSvs(RefGenomeInterface genome, String purpleDir, String linxDir, String refSample, String tumorSample,
            List<Sv> svs)
    {
        StringBuilder vcf = new StringBuilder();
        vcf.append("##fileformat=VCFv4.2\n");
        vcf.append("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"\">\n");
        vcf.append("##INFO=<ID=MATEID,Number=1,Type=String,Description=\"\">\n");
        vcf.append("##INFO=<ID=SVID,Number=1,Type=String,Description=\"\">\n");
        vcf.append("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"\">\n");
        vcf.append("##INFO=<ID=PURPLE_AF,Number=.,Type=Float,Description=\"\">\n");
        vcf.append("##INFO=<ID=PURPLE_CN,Number=.,Type=Float,Description=\"\">\n");
        vcf.append("##INFO=<ID=PURPLE_CN_CHANGE,Number=.,Type=Float,Description=\"\">\n");
        vcf.append("##INFO=<ID=PURPLE_JCN,Number=1,Type=Float,Description=\"\">\n");
        vcf.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">\n");
        vcf.append("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"\">\n");
        vcf.append("##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"\">\n");
        vcf.append("##FORMAT=<ID=VF,Number=1,Type=Integer,Description=\"\">\n");
        vcf.append("##FORMAT=<ID=REF,Number=1,Type=Integer,Description=\"\">\n");
        vcf.append("##FORMAT=<ID=REFPAIR,Number=1,Type=Integer,Description=\"\">\n");
        List<String> chromosomes = new ArrayList<>();
        svs.forEach(sv ->
        {
            chromosomes.add(sv.startChromosome());
            chromosomes.add(sv.endChromosome());
        });
        appendContigs(genome, vcf, chromosomes.stream().distinct().toList());
        vcf.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t").append(refSample).append('\t').append(tumorSample)
                .append('\n');

        List<LinxSvAnnotation> annotations = new ArrayList<>();
        List<LinxBreakend> breakends = new ArrayList<>();
        List<LinxFusion> fusions = new ArrayList<>();
        List<LinxDriver> drivers = new ArrayList<>();
        List<LinxCluster> clusters = new ArrayList<>();
        List<GeneCopyNumber> geneCopyNumbers = new ArrayList<>();

        int recordId = 100;
        int svId = 1;
        int breakendId = 1;
        int clusterId = 1;

        for(Sv sv : svs)
        {
            String startId = String.valueOf(recordId++);
            String endId = String.valueOf(recordId++);
            String svTag = startId + "_" + endId;

            String startRef = genome.getBaseString(sv.startChromosome(), sv.startPosition(), sv.startPosition()).toUpperCase();
            String endRef = genome.getBaseString(sv.endChromosome(), sv.endPosition(), sv.endPosition()).toUpperCase();

            String startAlt = bndAlt(
                    sv.startOrientation(), sv.endOrientation(), startRef, sv.endChromosome(), sv.endPosition(),
                    sv.insertSequence());
            String endAlt = bndAlt(sv.endOrientation(), sv.startOrientation(), endRef, sv.startChromosome(), sv.startPosition(), "");

            String svType = svTypeString(sv);

            appendSvRecord(
                    vcf, sv.startChromosome(), sv.startPosition(), startId, startRef, startAlt, svType, endId, svTag,
                    sv.purpleAf(), sv.copyNumberChange(), sv.junctionCopyNumber(), sv.tumorFragments());
            appendSvRecord(
                    vcf, sv.endChromosome(), sv.endPosition(), endId, endRef, endAlt, svType, startId, svTag,
                    sv.purpleAf(), sv.copyNumberChange(), sv.junctionCopyNumber(), sv.tumorFragments());

            int thisSvId = svId++;
            int thisClusterId = clusterId++;
            int startBreakendId = breakendId++;
            int endBreakendId = breakendId++;

            annotations.add(LinxTestFactory.svAnnotationBuilder()
                    .vcfIdStart(startId)
                    .vcfIdEnd(endId)
                    .svId(thisSvId)
                    .clusterId(thisClusterId)
                    .type(StructuralVariantType.valueOf(svType))
                    .build());

            ReportedStatus breakendStatus = sv.driverKind() == Sv.DriverKind.DISRUPTION
                    ? ReportedStatus.REPORTED
                    : ReportedStatus.NONE;
            breakends.add(LinxTestFactory.breakendBuilder()
                    .id(startBreakendId).svId(thisSvId).gene(sv.gene()).isStart(true).reportedStatus(breakendStatus).build());
            breakends.add(LinxTestFactory.breakendBuilder()
                    .id(endBreakendId).svId(thisSvId).gene(sv.gene()).isStart(false).reportedStatus(breakendStatus).build());

            clusters.add(ImmutableLinxCluster.builder()
                    .clusterId(thisClusterId).category("SIMPLE_SV").synthetic(false).resolvedType("SIMPLE_SV").clusterCount(1)
                    .clusterDesc("SV").build());

            switch(sv.driverKind())
            {
                case FUSION:
                    fusions.add(LinxTestFactory.fusionBuilder()
                            .name(sv.gene()).reported(true)
                            .fivePrimeBreakendId(startBreakendId).threePrimeBreakendId(endBreakendId).build());
                    break;
                case AMP:
                    drivers.add(LinxTestFactory.driverEventBuilder()
                            .clusterId(thisClusterId).eventType(DriverEventType.GAIN).gene(sv.gene()).build());
                    break;
                case DEL:
                    drivers.add(LinxTestFactory.driverEventBuilder()
                            .clusterId(thisClusterId).eventType(DriverEventType.DEL).gene(sv.gene()).build());
                    // Gene copy number whose minimum region is flanked by this SV's breakends (matches the DEL driver criteria).
                    geneCopyNumbers.add(geneCopyNumber(sv.gene(), sv.startChromosome(), sv.startPosition(), sv.endPosition()));
                    break;
                case DISRUPTION:
                    break;
            }
        }

        writeBgzip(PurpleCommon.purpleSomaticSvFile(purpleDir, tumorSample), vcf.toString());

        try
        {
            LinxSvAnnotation.write(LinxSvAnnotation.generateFilename(linxDir, tumorSample), annotations);
            LinxBreakend.write(LinxBreakend.generateFilename(linxDir, tumorSample), breakends);
            LinxFusion.write(LinxFusion.generateFilename(linxDir, tumorSample), fusions);
            LinxDriver.write(LinxDriver.generateFilename(linxDir, tumorSample), drivers);
            LinxCluster.write(LinxCluster.generateFilename(linxDir, tumorSample), clusters);
            GeneCopyNumberFile.write(GeneCopyNumberFile.generateFilename(purpleDir, tumorSample), geneCopyNumbers);
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    private static void appendSvRecord(StringBuilder vcf, String chromosome, int position, String id, String ref, String alt,
            String svType, String mateId, String svTag, double af, double cnChange, double jcn, int tumorFragments)
    {
        String info = format(
                "SVTYPE=%s;MATEID=%s;SVID=%s;CIPOS=0,0;PURPLE_AF=%.3f,%.3f;PURPLE_CN=2.0,2.0;PURPLE_CN_CHANGE=%.3f,%.3f;PURPLE_JCN=%.3f",
                svType, mateId, svTag, af, af, cnChange, cnChange, jcn);
        String refGenotype = "0/0:30,0:0.0:0:30:10";
        String tumorGenotype = format("0/1:40,%d:%.3f:%d:30:15", tumorFragments, af, tumorFragments);
        vcf.append(chromosome).append('\t')
                .append(position).append('\t')
                .append(id).append('\t')
                .append(ref).append('\t')
                .append(alt).append('\t')
                .append("100").append('\t')
                .append("PASS").append('\t')
                .append(info).append('\t')
                .append("GT:AD:AF:VF:REF:REFPAIR").append('\t')
                .append(refGenotype).append('\t')
                .append(tumorGenotype).append('\n');
    }

    private static void writeGermlineSvs(RefGenomeInterface genome, String linxGermlineDir, String sampleId, List<GermlineSv> svs)
    {
        List<LinxGermlineDisruption> disruptions = new ArrayList<>();
        List<LinxBreakend> breakends = new ArrayList<>();

        int svId = 1;
        int breakendId = 1;
        for(GermlineSv sv : svs)
        {
            int thisSvId = svId++;
            disruptions.add(new LinxGermlineDisruption(
                    thisSvId, String.valueOf(thisSvId), sv.startChromosome(), sv.endChromosome(),
                    sv.startPosition(), sv.endPosition(), sv.startOrientation(), sv.endOrientation(), StructuralVariantType.BND,
                    "PASS", "", 100.0,
                    "", "", 0, 0,
                    1.0,
                    0.3, 0.3, 2.0, 2.0,
                    1.0, 1.0,
                    20, 10, 10,
                    20, 10, 10,
                    sv.insertSequence(), "", "", "",
                    sv.gene(), 1, 1, "SIMPLE_SV",
                    "", "", 0));

            breakends.add(LinxTestFactory.breakendBuilder()
                    .id(breakendId++).svId(thisSvId).gene(sv.gene()).isStart(true).reportedStatus(ReportedStatus.REPORTED).build());
            breakends.add(LinxTestFactory.breakendBuilder()
                    .id(breakendId++).svId(thisSvId).gene(sv.gene()).isStart(false).reportedStatus(ReportedStatus.REPORTED).build());
        }

        try
        {
            LinxGermlineDisruption.write(LinxGermlineDisruption.generateFilename(linxGermlineDir, sampleId), disruptions);
            LinxBreakend.write(LinxBreakend.generateFilename(linxGermlineDir, sampleId, true), breakends);
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    private static String bndAlt(byte thisOrientation, byte otherOrientation, String refBase, String mateChromosome, int matePosition,
            String insert)
    {
        String bracket = otherOrientation == ORIENT_FWD ? "]" : "[";
        String locus = bracket + mateChromosome + ":" + matePosition + bracket;
        if(thisOrientation == ORIENT_FWD)
        {
            return refBase + insert + locus;
        }
        else
        {
            return locus + insert + refBase;
        }
    }

    private static String svTypeString(Sv sv)
    {
        if(!sv.startChromosome().equals(sv.endChromosome()))
        {
            return "BND";
        }
        if(sv.startOrientation() == ORIENT_FWD && sv.endOrientation() == ORIENT_REV)
        {
            return "DEL";
        }
        if(sv.startOrientation() == ORIENT_REV && sv.endOrientation() == ORIENT_FWD)
        {
            return "DUP";
        }
        return "INV";
    }

    private static GeneCopyNumber geneCopyNumber(String gene, String chromosome, int minRegionStart, int minRegionEnd)
    {
        return new GeneCopyNumber(
                chromosome, minRegionStart, minRegionEnd, gene, "", true, "p1", 2.0, 0.5,
                0.0, 1, 1, minRegionStart, minRegionEnd, 1, 0.45,
                SegmentSupport.NONE, SegmentSupport.NONE, CopyNumberMethod.UNKNOWN, 0.25);
    }

    private static String impactString(String gene, CodingEffect codingEffect)
    {
        VariantImpact impact = new VariantImpact(
                gene, "ENST00000000001", codingEffect == CodingEffect.NONE ? "intron_variant" : "missense_variant",
                codingEffect, "", "", false, "", codingEffect, gene.isEmpty() ? 0 : 1);
        return String.join(",", VariantImpactSerialiser.toVcfData(impact));
    }

    private static void appendContigs(RefGenomeInterface genome, StringBuilder vcf, List<String> chromosomes)
    {
        for(String chromosome : chromosomes)
        {
            int length = genome.getChromosomeLength(chromosome);
            vcf.append(format("##contig=<ID=%s,length=%d>\n", chromosome, length));
        }
    }

    private static String differentBase(String base)
    {
        return switch(base.toUpperCase())
        {
            case "A" -> "T";
            case "T" -> "A";
            case "C" -> "G";
            case "G" -> "C";
            default -> "A";
        };
    }

    private static void writeBgzip(String path, String content)
    {
        try(BlockCompressedOutputStream out = new BlockCompressedOutputStream(new File(path)))
        {
            out.write(content.getBytes(StandardCharsets.US_ASCII));
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }

    private static void createDir(String path)
    {
        try
        {
            Files.createDirectories(Path.of(path));
        }
        catch(IOException e)
        {
            throw new UncheckedIOException(e);
        }
    }
}
