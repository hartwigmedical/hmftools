package com.hartwig.hmftools.isofox.cram;

import com.hartwig.computeengine.execution.vm.VmDirectories;

import static java.lang.String.format;

public enum ExternalTool {
    BAMCOMP("bamcomp", "bamcomp.jar", "1.3"),
    SAMBAMBA("sambamba", "sambamba", "0.6.8"),
    SAMTOOLS("samtools", "samtools", "1.14");

    private final String toolName;
    private final String version;
    private final String binary;

    ExternalTool(final String toolName, final String binary, final String version) {
        this.toolName = toolName;
        this.version = version;
        this.binary = binary;
    }

    public String path() {
        return format("%s/%s/%s", VmDirectories.TOOLS, toolName, version);
    }

    public String binaryPath() {
        return format("%s/%s/%s/%s", VmDirectories.TOOLS, toolName, version, binary);
    }

    public String getToolName() {
        return toolName;
    }

    public String getVersion() {
        return version;
    }

    public String getBinary() {
        return binary;
    }
}
