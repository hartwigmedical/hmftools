package com.hartwig.hmftools.isofox.cram;



import java.util.List;

import com.google.common.collect.ImmutableList;
import com.hartwig.computeengine.execution.vm.Bash;
import com.hartwig.computeengine.execution.vm.command.BashCommand;

import static com.hartwig.hmftools.isofox.cram.ExternalTool.BAMCOMP;
import static com.hartwig.hmftools.isofox.cram.ExternalTool.SAMBAMBA;
import static com.hartwig.hmftools.isofox.cram.ExternalTool.SAMTOOLS;

public class CramAndValidateCommands {
    String RESOURCES = "/opt/resources";
    String REFERENCE_GENOME = "reference_genome";
    String REF_FILE = "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna";

    private final String inputBam;
    private final String outputCram;



    public CramAndValidateCommands(final String inputBam, final String outputCram) {
        this.inputBam = inputBam;
        this.outputCram = outputCram;
    }

    public List<BashCommand> commands() {
        return ImmutableList.of(new VersionedToolCommand(SAMTOOLS.getToolName(),
                        SAMTOOLS.getBinary(),
                        SAMTOOLS.getVersion(),
                        "view",
                        "-T",
                        refGenomeFile(),
                        "-o",
                        outputCram,
                        "-O",
                        "cram,embed_ref=1",
                        "-@",
                        Bash.allCpus(),
                        inputBam),
                new VersionedToolCommand(SAMTOOLS.getToolName(),
                        SAMTOOLS.getBinary(),
                        SAMTOOLS.getVersion(),
                        "reheader",
                        "--no-PG",
                        "--in-place",
                        "--command",
                        "'grep -v ^@PG'",
                        outputCram),
                new VersionedToolCommand(SAMTOOLS.getToolName(), SAMTOOLS.getBinary(), SAMTOOLS.getVersion(), "index", outputCram),
                JavaCommandFactory.javaClassCommand(BAMCOMP, "com.hartwig.bamcomp.BamCompMain", "4G", javaClassCommandArguments()));
    }

    private List<String> javaClassCommandArguments() {
        return List.of("-r",
                refGenomeFile(),
                "-1",
                inputBam,
                "-2",
                outputCram,
                "-n",
                "6",
                "--samtools-binary",
                SAMTOOLS.binaryPath(),
                "--sambamba-binary",
                SAMBAMBA.binaryPath());
    }

    private String refGenomeFile() {
        return String.format("%s/%s/%s/%s", RESOURCES, REFERENCE_GENOME, "37", REF_FILE);
    }

}
