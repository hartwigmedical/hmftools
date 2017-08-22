package com.hartwig.hmftools.civic;

import com.hartwig.hmftools.civic.api.CivicApiWrapper;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class CivicClientApplication {

    private static final Logger LOGGER = LogManager.getLogger(CivicClientApplication.class);

    private static final String SNP_VCF_LINE =
            "7\t55249071\trs121434569;COSM6240\tC\tT\t161.76\tPASS\tAB=0.123288;ABP=272.957;AC=3;AF=0.273;AN=11;ANN=T|missense_variant|MODERATE|EGFR|ENSG00000146648|transcript|ENST00000275493|protein_coding|20/28|c.2369C>T|p.Thr790Met|2546/9821|2369/3633|790/1210||,T|missense_variant|MODERATE|EGFR|ENSG00000146648|transcript|ENST00000455089|protein_coding|19/26|c.2234C>T|p.Thr745Met|2491/3844|2234/3276|745/1091||,T|missense_variant|MODERATE|EGFR|ENSG00000146648|transcript|ENST00000454757|protein_coding|20/28|c.2210C>T|p.Thr737Met|2393/3938|2210/3474|737/1157||,T|sequence_feature|LOW|EGFR|ENSG00000146648|topological_domain:Cytoplasmic|ENST00000275493|protein_coding||c.2369C>T||||||,T|sequence_feature|LOW|EGFR|ENSG00000146648|domain:Protein_kinase|ENST00000275493|protein_coding||c.2369C>T||||||,T|sequence_feature|LOW|EGFR|ENSG00000146648|beta_strand|ENST00000275493|protein_coding||c.2369C>T||||||,T|sequence_feature|LOW|EGFR|ENSG00000146648|nucleotide_phosphate-binding_region:ATP|ENST00000275493|protein_coding||c.2369C>T||||||,T|intron_variant|MODIFIER|EGFR|ENSG00000146648|transcript|ENST00000442591|protein_coding|17/17|c.*28+8450C>T||||||,T|non_coding_exon_variant|MODIFIER|EGFR-AS1|ENSG00000224057|transcript|ENST00000442411|antisense|2/2|n.1178G>A||||||;AO=27;CIGAR=1X;DB;DP=244;DPB=244.0;DPRA=8.76;EPP=3.09072;EPPR=3.0103;GTI=0;LEN=1;MEANALT=2.0;MQM=60.0;MQMR=60.0;NS=2;NT=ref;NUMALT=1;ODDS=18.3095;PAIRED=1.0;PAIREDR=1.0;PAO=0.0;PQA=0.0;PQR=0.0;PRO=0.0;QA=950;QR=7901;QSS=14;QSS_NT=14;RO=216;RPL=11.0;RPP=5.02092;RPPR=6.26751;RPR=16.0;RUN=1;SAF=15;SAP=3.73412;SAR=12;SGT=CC->CT;SOMATIC;SRF=93;SRP=12.0581;SRR=123;TQSS=1;TQSS_NT=1;TYPE=snp;VT=SNP;set=freebayes-mutect-strelka;technology.ILLUMINA=1.0;CSA=3,3;CSP=3\tGT:AD:DP\t0/1:180,24:208";
    private static final String INDEL_VCF_LINE =
            "7\t55242464\trs121913421;COSM6223\tAGGAATTAAGAGAAGC\tA\t1656.44\tPASS\tAB=0.472527;ABP=4.20342;AC=2;AF=0.25;AN=8;ANN=A|disruptive_inframe_deletion|MODERATE|EGFR|ENSG00000146648|transcript|ENST00000275493|protein_coding|19/28|c.2235_2249delGGAATTAAGAGAAGC|p.Glu746_Ala750del|2412/9821|2235/3633|745/1210||,A|disruptive_inframe_deletion|MODERATE|EGFR|ENSG00000146648|transcript|ENST00000455089|protein_coding|18/26|c.2100_2114delGGAATTAAGAGAAGC|p.Glu701_Ala705del|2357/3844|2100/3276|700/1091||,A|disruptive_inframe_deletion|MODERATE|EGFR|ENSG00000146648|transcript|ENST00000454757|protein_coding|19/28|c.2076_2090delGGAATTAAGAGAAGC|p.Glu693_Ala697del|2259/3938|2076/3474|692/1157||,A|sequence_feature|LOW|EGFR|ENSG00000146648|topological_domain:Cytoplasmic|ENST00000275493|protein_coding||c.2235_2249delGGAATTAAGAGAAGC||||||,A|sequence_feature|LOW|EGFR|ENSG00000146648|domain:Protein_kinase|ENST00000275493|protein_coding||c.2235_2249delGGAATTAAGAGAAGC||||||,A|sequence_feature|LOW|EGFR|ENSG00000146648|beta_strand|ENST00000275493|protein_coding||c.2235_2249delGGAATTAAGAGAAGC||||||,A|sequence_feature|LOW|EGFR|ENSG00000146648|binding_site:ATP|ENST00000275493|protein_coding||c.2235_2249delGGAATTAAGAGAAGC||||||,A|intron_variant|MODIFIER|EGFR|ENSG00000146648|transcript|ENST00000442591|protein_coding|17/17|c.*28+1844_*28+1858delGGAATTAAGAGAAGC||||||;AO=86;CIGAR=1M15D1M;DB;DP=208;DPB=132.118;DPRA=7.0;EPP=20.0791;EPPR=6.20142;GTI=0;IC=0;IHP=3;LEN=15;MEANALT=3.0;MQM=60.1163;MQMR=60.0;NS=2;NT=ref;NUMALT=1;ODDS=17.6164;OLD_VARIANT=7:55242464:AGGAATTAAGAGAAGCA/AA;PAIRED=0.988372;PAIREDR=1.0;PAO=0.0;PQA=0.0;PQR=0.0;PRO=0.0;QA=2531;QR=4416;QSI=65;QSI_NT=65;RC=1;RO=115;RPL=42.0;RPP=3.1113;RPPR=3.02918;RPR=44.0;RU=GGAATTAAGAGAAGC;RUN=1;SAF=60;SAP=32.1989;SAR=26;SGT=ref->het;SOMATIC;SRF=53;SRP=4.53977;SRR=62;TQSI=1;TQSI_NT=1;TYPE=del;set=freebayes-strelka;technology.ILLUMINA=1.0;CSA=2,2;CSP=2\tGT:AD:DP\t0/1:106,85:196";

    public static void main(final String... args) throws Exception {
        final Options options = createOptions();
        final CommandLine cmd = createCommandLine(options, args);

        final SomaticVariant snpVariant = SomaticVariantFactory.fromVCFLine(SNP_VCF_LINE);
        final SomaticVariant indelVariant = SomaticVariantFactory.fromVCFLine(INDEL_VCF_LINE);

        CivicApiWrapper.getVariantsAtPosition(1956, snpVariant)
                //                .flatMapIterable(CivicVariant::evidenceItems)
                //                .filter(evidenceItem -> !evidenceItem.drugs().isEmpty())
                //                .filter(evidenceItem -> evidenceItem.level() <= 'C')
                .subscribe(LOGGER::info, Throwable::printStackTrace, () -> LOGGER.info("completed"));
    }

    @NotNull
    private static Options createOptions() {
        final Options options = new Options();
        return options;
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final Options options, @NotNull final String... args) throws ParseException {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
