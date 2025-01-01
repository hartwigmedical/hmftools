package com.hartwig.hmftools.isofox.cram;

import java.util.Collections;
import java.util.List;

import com.hartwig.computeengine.execution.vm.command.java.JavaClassCommand;
import com.hartwig.computeengine.execution.vm.command.java.JavaJarCommand;

public final class JavaCommandFactory {

    private JavaCommandFactory() {
    }

    public static JavaJarCommand javaJarCommand(HmfTool hmfTool, List<String> arguments) {
        return new JavaJarCommand(hmfTool.getToolName(), hmfTool.runVersion(), hmfTool.jar(), hmfTool.maxHeapStr(), arguments);
    }

    public static JavaClassCommand javaClassCommand(HmfTool hmfTool, String mainClass, List<String> arguments) {
        return new JavaClassCommand(hmfTool.getToolName(),
                hmfTool.runVersion(),
                hmfTool.jar(),
                mainClass,
                hmfTool.maxHeapStr(),
                Collections.emptyList(),
                arguments);
    }

    public static JavaClassCommand javaClassCommand(ExternalTool externalTool, String mainClass, String heapSize, List<String> arguments) {
        return new JavaClassCommand(externalTool.getToolName(),
                externalTool.getVersion(),
                externalTool.getBinary(),
                mainClass,
                heapSize,
                Collections.emptyList(),
                arguments);
    }
}
