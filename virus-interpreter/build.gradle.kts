
plugins {
    id("com.hartwig.java-conventions")
}

description = "HMF Tools - Virus Interpreter"

dependencies {
    implementation(project(":hmf-common"))
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.junit)
    testImplementation(project(path = ":hmf-common", configuration = "tests"))
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.virusinterpreter.VirusInterpreterApplication")
}
