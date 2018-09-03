package com.hartwig.hmftools.portal.converter

import com.google.common.io.Resources
import com.hartwig.hmftools.common.context.ProductionRunContextFactory
import io.kotlintest.matchers.shouldBe
import io.kotlintest.specs.StringSpec
import java.io.File
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.util.stream.Collectors

class PortalDataConverterApplicationKtTest : StringSpec() {
    private val patientsDataFilePath = Resources.getResource("patient_data_file.csv").path
    private val runsV3 = readRun(Paths.get(Resources.getResource("runs_v3.txt").path))
    private val runsV4 = readRun(Paths.get(Resources.getResource("runs_v4.txt").path))

    init {
        "can import pipeline v3 vcf"{
            val testDir = createTestDir()
            convertSamples(runsV3, testDir.toString(), File(patientsDataFilePath))
            outputAsExpected(testDir) shouldBe true
            deleteTestDir(testDir)
        }

        "can import pipeline v4 vcf"{
            val testDir = createTestDir()
            convertSamples(runsV4, testDir.toString(), File(patientsDataFilePath))
            outputAsExpected(testDir) shouldBe true
            deleteTestDir(testDir)
        }
    }

    private fun outputAsExpected(path: Path): Boolean {
        val expectedOutputPath = Paths.get(Resources.getResource("expectedOutput").path)
        val expectedFilesAndDirs = subFilesAndDirs(expectedOutputPath)
        val filesAndDirs = subFilesAndDirs(path)
        return filesEqual(filesAndDirs, expectedFilesAndDirs)
    }

    private fun subFilesAndDirs(path: Path) = path.toFile().walkTopDown().filterNot { it.path == path.toString() }.toList()

    private fun filesEqual(files: List<File>, others: List<File>) =
            files.all { fileInFiles(it, others) } && others.all { fileInFiles(it, files) }

    private fun fileInFiles(file: File, files: List<File>): Boolean {
        val filesWithSameName = files.filter { it.name == file.name }
        return if (file.isFile) filesWithSameName.size == 1 && linesEqual(file, filesWithSameName[0])
        else filesWithSameName.size == 1
    }

    private fun linesEqual(file: File, other: File) = file.readLines() == other.readLines()

    private fun readRun(path: Path) = Files.lines(path).map { it.trim() }.filter { it.isNotEmpty() }.collect(Collectors.toList()).map {
        ProductionRunContextFactory.fromRunDirectory(Resources.getResource(it).path)
    }

    private fun createTestDir() = Files.createTempDirectory(Paths.get(Resources.getResource(".").path), "portal")

    private fun deleteTestDir(path: Path) = path.toFile().deleteRecursively()
}
