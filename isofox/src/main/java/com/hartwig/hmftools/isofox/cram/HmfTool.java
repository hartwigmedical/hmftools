package com.hartwig.hmftools.isofox.cram;


import static java.lang.String.format;

import com.hartwig.computeengine.execution.vm.VmDirectories;

public enum HmfTool {

    AMBER("4.0.1", 20, 24, 16, false),
    BAM_TOOLS("1.2.1", 16, 24, 16, false),
    CHORD("2.02_1.14", Defaults.JAVA_HEAP, 12, 4, false),
    CIDER("1.0.3", 16, 24, 4, false),
    COBALT("1.16", 20, 24, 16, false),
    CUPPA("2.3.0", Defaults.JAVA_HEAP, 16, 4, false),
    GRIDSS("2.13.3", Defaults.JAVA_HEAP, 64, 24, false),
    GRIPSS("2.4", 16, 24, 4, false),
    HEALTH_CHECKER("3.5", Defaults.JAVA_HEAP, 32, 8, false),
    LILAC("1.6", 16, 24, 8, false),
    LINX("1.25", 8, 12, 4, false),
    MARK_DUPS("1.1.7", 40, 120, 24, false),
    ORANGE("3.7.0", 16, 18, 4, false),
    PAVE("1.6", 30, 40, 8, false),
    PEACH("1.8", 1, 4, 2, false),
    PURPLE("4.0.2", 30, 40, 8, false),
    SAGE("3.4.3", 48, 64, 24, false),
    SIGS("1.2.1", Defaults.JAVA_HEAP, 16, 4, false),
    SV_PREP("1.2.4", 48, 64, 24, false),
    TEAL("1.2.2", 30, 32, 16, false),
    VIRUSBREAKEND_GRIDSS("2.13.3", Defaults.JAVA_HEAP, 64, 12, false),
    VIRUS_INTERPRETER("1.5.0", Defaults.JAVA_HEAP, 8, 2, false);

    private static final String PILOT_VERSION = "pilot"; // will pick up the jar from /opt/toolName/pilot/toolName.jar
    private final String toolName;
    private final String version;
    private final int maxHeap;
    private final int memoryGb;
    private final int cpus;
    private final boolean usePilot;

    HmfTool(final String version, final int maxHeap, final int memoryGb, final int cpus, final boolean usePilot) {
        toolName = this.toString().toLowerCase().replace('_', '-');
        this.version = version;
        this.maxHeap = maxHeap;
        this.memoryGb = memoryGb;
        this.cpus = cpus;
        this.usePilot = usePilot;
    }

    public String getToolName() {
        return toolName;
    }

    public String runVersion() {
        return usePilot ? PILOT_VERSION : version;
    }

    public String versionInfo() {
        return usePilot ? format("%s, using pilot", version) : version;
    }

    public int getMaxHeap() {
        return maxHeap;
    }

    public int getMemoryGb() {
        return memoryGb;
    }

    public int getCpus() {
        return cpus;
    }

    public String directory() {
        return toolName;
    }

    public String jar() {
        return format("%s.jar", toolName);
    }

    public String maxHeapStr() {
        return format("%dG", maxHeap);
    }

    public String jarPath() {
        return format("%s/%s/%s/%s", VmDirectories.TOOLS, directory(), version, jar());
    }

    private static final class Defaults {
        private static final int JAVA_HEAP = 4;
    }
}