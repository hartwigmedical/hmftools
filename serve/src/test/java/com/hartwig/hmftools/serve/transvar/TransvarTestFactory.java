package com.hartwig.hmftools.serve.transvar;

import java.io.FileNotFoundException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;

import org.jetbrains.annotations.NotNull;

public final class TransvarTestFactory {

    private static final String REF_GENOME_FASTA_FILE = Resources.getResource("refgenome/v37/ref.fasta").getPath();

    private TransvarTestFactory() {
    }

    @NotNull
    static Transvar testTransvar(@NotNull TransvarProcess process) {
        return new Transvar(process, testInterpreter(), HmfGenePanelSupplier.allGenesMap37());
    }

    @NotNull
    static TransvarInterpreter testInterpreter() {
        try {
            return TransvarInterpreter.fromRefGenomeFastaFile(REF_GENOME_FASTA_FILE);
        } catch (FileNotFoundException exception) {
            throw new IllegalStateException("Cannot create test interpreter! Message=" + exception.getMessage());
        }
    }
}
