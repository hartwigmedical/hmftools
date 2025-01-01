package com.hartwig.hmftools.isofox.cram;

import java.util.Arrays;
import java.util.List;

import com.hartwig.computeengine.execution.vm.VmDirectories;
import com.hartwig.computeengine.execution.vm.command.BashCommand;

public class  VersionedToolCommand implements BashCommand {

    private final String toolName;
    private final String toolBinaryName;
    private final String version;
    private final List<String> arguments;

    public VersionedToolCommand(final String toolName, final String toolBinaryName, final String version, final List<String> arguments) {
        this.toolName = toolName;
        this.toolBinaryName = toolBinaryName;
        this.version = version;
        this.arguments = arguments;
    }

    public VersionedToolCommand(final String toolName, final String toolBinaryName, final String version, final String... arguments) {
        this(toolName, toolBinaryName, version, Arrays.asList(arguments));
    }

    @Override
    public String asBash() {
        return String.format("%s/%s/%s/%s %s",
                VmDirectories.TOOLS,
                toolName,
                version,
                toolBinaryName,
                String.join(" ", arguments));
    }
}