> **Example Note:** This is a whole-codebase example for demonstration purposes.
> In typical usage, arch lens diagrams are scoped to the subsystem being
> modified/added/removed by a plan — not the entire project. The result is
> normally much simpler and more focused than what you see here.

# Development Diagram: STAR Aligner

**Lens:** Development (Development)
**Question:** How is it built and tested?
**Date:** 2026-02-14
**Scope:** Full STAR aligner build and development workflow

## Overview

| Aspect | Description |
|--------|-------------|
| **Build System** | GNU Make with platform-specific targets |
| **Language** | C++11 with some C (zlib, htslib interfaces) |
| **Compiler Support** | GCC 4.8+, Clang 3.4+, Intel C++ |
| **Optimization** | -O3, -march=native, optional AVX2/SSE4.2 |
| **Test Strategy** | Regression tests with known datasets, output comparison |
| **Source Organization** | Flat source/ directory with module prefixes |

## Development Workflow

```mermaid
%%{init: {'flowchart': {'nodeSpacing': 50, 'rankSpacing': 60, 'curve': 'basis'}}}%%
flowchart TB
    Start["Developer Starts Build<br/>━━━━━━━━━━<br/>Clone STAR repo<br/>git clone STAR"]:::terminal

    subgraph SourceCode["Source Code Organization"]
        SourceDir["source/ Directory<br/>━━━━━━━━━━<br/>~150 .cpp/.h files<br/>Flat structure"]:::cli

        subgraph Modules["Code Modules (by prefix)"]
            Genome["Genome*.cpp/h<br/>━━━━━━━━━━<br/>Index loading, management<br/>GenomeGenerate.cpp"]:::newComponent

            Align["align*.cpp/h<br/>━━━━━━━━━━<br/>Core alignment logic<br/>alignReads.cpp"]:::newComponent

            SAM["BAM*.cpp/h, SAM*.cpp/h<br/>━━━━━━━━━━<br/>Output formatting<br/>BAMoutput.cpp"]:::newComponent

            Params["Parameters*.cpp/h<br/>━━━━━━━━━━<br/>Command-line parsing<br/>ParametersChimeric.cpp"]:::newComponent
        end

        Headers["Header Files<br/>━━━━━━━━━━<br/>Class declarations<br/>Inline functions"]:::stateNode
    end

    subgraph BuildConfig["Build Configuration"]
        Makefile["Makefile<br/>━━━━━━━━━━<br/>Top-level build rules<br/>Platform detection"]:::handler

        PlatformTargets{"Platform<br/>Target"}:::phase

        LinuxTarget["Linux Build Config<br/>━━━━━━━━━━<br/>-lpthread -lz -lrt<br/>GCC defaults"]:::newComponent

        MacTarget["macOS Build Config<br/>━━━━━━━━━━<br/>Clang defaults<br/>Different linker flags"]:::newComponent

        StaticTarget["Static Build Config<br/>━━━━━━━━━━<br/>-static linking<br/>Portable binary"]:::newComponent
    end

    subgraph Compilation["Compilation Process"]
        CXXFlags["Compiler Flags<br/>━━━━━━━━━━<br/>-O3 -std=c++11<br/>-pthread -march=native"]:::handler

        OptionalFlags{"Optional<br/>Optimizations"}:::phase

        AVX2["AVX2 SIMD<br/>━━━━━━━━━━<br/>-mavx2 -mfma<br/>Faster alignment"]:::newComponent

        SSE4["SSE4.2 SIMD<br/>━━━━━━━━━━<br/>-msse4.2<br/>Baseline acceleration"]:::newComponent

        Htslib["HTSlib Support<br/>━━━━━━━━━━<br/>-DUSE_HTSLIB<br/>Link against htslib"]:::integration

        CompileUnits["Compile .cpp Files<br/>━━━━━━━━━━<br/>g++ -c source/*.cpp<br/>Generate .o files"]:::handler

        LinkBinary["Link STAR Binary<br/>━━━━━━━━━━<br/>g++ *.o -lz -lpthread<br/>Produce executable"]:::handler
    end

    subgraph Testing["Testing Workflow"]
        TestData["Test Datasets<br/>━━━━━━━━━━<br/>Small genome + reads<br/>Known outputs"]:::cli

        RunIndex["Generate Test Index<br/>━━━━━━━━━━<br/>./STAR --runMode genomeGenerate<br/>Validate index files"]:::output

        RunAlign["Run Test Alignment<br/>━━━━━━━━━━<br/>./STAR --runMode alignReads<br/>Align test reads"]:::output

        CompareOutput["Compare Outputs<br/>━━━━━━━━━━<br/>diff expected vs actual<br/>SAM/BAM validation"]:::detector

        TestPass{"All Tests<br/>Pass?"}:::phase
    end

    subgraph DebugBuild["Debug Build (Optional)"]
        DebugFlags["Debug Compiler Flags<br/>━━━━━━━━━━<br/>-g -O0 -DDEBUG<br/>No optimization"]:::handler

        DebugBinary["STAR Debug Binary<br/>━━━━━━━━━━<br/>Large binary with symbols<br/>gdb/lldb compatible"]:::output
    end

    subgraph ReleaseArtifacts["Release Artifacts"]
        ReleaseBinary["STAR Binary<br/>━━━━━━━━━━<br/>Optimized executable<br/>~2-3 MB stripped"]:::output

        Documentation["Documentation<br/>━━━━━━━━━━<br/>doc/STARmanual.pdf<br/>Parameter descriptions"]:::terminal

        SourceRelease["Source Tarball<br/>━━━━━━━━━━<br/>STAR_x.y.z.tar.gz<br/>GitHub release"]:::terminal
    end

    Start --> SourceDir
    SourceDir --> Genome
    SourceDir --> Align
    SourceDir --> SAM
    SourceDir --> Params
    SourceDir --> Headers

    SourceDir --> Makefile
    Makefile --> PlatformTargets

    PlatformTargets -->|"make STAR"| LinuxTarget
    PlatformTargets -->|"make STARforMacStatic"| MacTarget
    PlatformTargets -->|"make STARstatic"| StaticTarget

    LinuxTarget --> CXXFlags
    MacTarget --> CXXFlags
    StaticTarget --> CXXFlags

    CXXFlags --> OptionalFlags

    OptionalFlags -->|"CXXFLAGS+=-mavx2"| AVX2
    OptionalFlags -->|"Default SSE4.2"| SSE4
    OptionalFlags -->|"CXXFLAGS+=-DUSE_HTSLIB"| Htslib
    OptionalFlags -->|"Standard build"| CompileUnits

    AVX2 --> CompileUnits
    SSE4 --> CompileUnits
    Htslib --> CompileUnits

    CompileUnits --> LinkBinary

    LinkBinary --> TestData
    TestData --> RunIndex
    RunIndex --> RunAlign
    RunAlign --> CompareOutput
    CompareOutput --> TestPass

    TestPass -->|"No"| DebugFlags
    TestPass -->|"Yes"| ReleaseBinary

    DebugFlags --> DebugBinary
    DebugBinary --> CompareOutput

    ReleaseBinary --> Documentation
    ReleaseBinary --> SourceRelease

    classDef cli fill:#1a237e,stroke:#7986cb,stroke-width:2px,color:#fff;
    classDef stateNode fill:#004d40,stroke:#4db6ac,stroke-width:2px,color:#fff;
    classDef handler fill:#e65100,stroke:#ffb74d,stroke-width:2px,color:#fff;
    classDef phase fill:#6a1b9a,stroke:#ba68c8,stroke-width:2px,color:#fff;
    classDef newComponent fill:#2e7d32,stroke:#81c784,stroke-width:2px,color:#fff;
    classDef output fill:#00695c,stroke:#4db6ac,stroke-width:2px,color:#fff;
    classDef detector fill:#b71c1c,stroke:#ef5350,stroke-width:2px,color:#fff;
    classDef gap fill:#ff6f00,stroke:#ffa726,stroke-width:2px,color:#000;
    classDef integration fill:#c62828,stroke:#ef9a9a,stroke-width:2px,color:#fff;
    classDef terminal fill:#1a237e,stroke:#7986cb,stroke-width:2px,color:#fff;
```

## Color Legend

- **CLI (Dark Blue)**: Source code entry points and test data
- **Handler (Orange)**: Build system and compiler configuration
- **Phase (Purple)**: Build decision points
- **New Component (Green)**: Source code modules and build variants
- **State Node (Teal)**: Code organization structures
- **Output (Dark Teal)**: Build outputs and test results
- **Detector (Red)**: Validation and comparison steps
- **Integration (Dark Red)**: External library integration
- **Terminal (Dark Blue)**: Final release artifacts

## Analysis

### Source Code Organization

STAR uses a **flat source directory** with module prefixes:

```
source/
├── STAR.cpp                    # Main entry point
├── Genome*.cpp/h               # Genome index (15 files)
├── Parameters*.cpp/h           # Command-line parsing (20 files)
├── alignReads*.cpp/h           # Alignment engine (30 files)
├── BAMoutput*.cpp/h            # BAM writing (10 files)
├── SAMoutput*.cpp/h            # SAM writing (5 files)
├── ReadAlign*.cpp/h            # Read alignment classes (25 files)
├── SuffixArray*.cpp/h          # Suffix array construction (8 files)
├── Transcript*.cpp/h           # Transcript handling (12 files)
└── [other modules]
```

**Why flat structure?**
- Simpler build system (all files at same level)
- Easier to search and navigate
- Avoids deep include paths
- Module prefixes provide logical grouping

### Build System Architecture

The Makefile uses **target-based builds** for different platforms:

```makefile
# Core variables
CXX := g++
CXXFLAGS := -O3 -std=c++11 -pthread -march=native
LDFLAGS := -lz -lrt

# Platform-specific targets
STAR: # Linux dynamic
    $(CXX) $(CXXFLAGS) source/*.cpp -o STAR $(LDFLAGS)

STARstatic: # Linux static
    $(CXX) $(CXXFLAGS) -static source/*.cpp -o STAR $(LDFLAGS)

STARforMacStatic: # macOS
    clang++ $(CXXFLAGS) -Xpreprocessor -fopenmp source/*.cpp -o STAR -lz
```

### Compiler Flags Breakdown

**Standard flags** (always used):
- `-O3`: Aggressive optimization (loop unrolling, inlining)
- `-std=c++11`: C++11 standard (range-based for, auto, lambdas)
- `-pthread`: POSIX threads support
- `-march=native`: Optimize for build machine CPU

**Optional flags** (user can add):
- `-mavx2 -mfma`: AVX2 SIMD for vectorized alignment
- `-DUSE_HTSLIB`: Enable htslib for BAM output
- `-g -O0`: Debug build with symbols
- `-static`: Static linking for portability

### Platform Differences

**Linux**:
- Uses GCC by default
- Links against `-lrt` (realtime library for shm_open)
- Requires `-lpthread` explicitly

**macOS**:
- Uses Clang (GCC compatibility mode)
- No `-lrt` needed (shm_open in libc)
- OpenMP requires `-Xpreprocessor -fopenmp` and separate libomp
- Static linking more complex (macOS discourages it)

**Windows**:
- Not officially supported
- Community builds via Cygwin or WSL
- Some success with MinGW-w64

### SIMD Optimizations

STAR has hand-optimized SIMD code for alignment:

1. **AVX2** (256-bit vectors):
   - Process 32 bases in parallel
   - Smith-Waterman scoring
   - Requires `-mavx2 -mfma`
   - ~30% faster on Haswell+ CPUs

2. **SSE4.2** (128-bit vectors):
   - Baseline SIMD support
   - Compatible with most x86-64 CPUs
   - Used if AVX2 not available

3. **Auto-detection**:
   - Runtime CPU feature detection
   - Falls back to scalar code if SIMD unavailable

### Dependency Management

STAR has minimal external dependencies:

**Required**:
- **zlib**: Included in most systems, statically linkable
- **pthreads**: Part of libc on POSIX systems
- **OpenMP**: libgomp (GCC) or libomp (Clang)

**Optional**:
- **htslib**: For professional BAM/CRAM output
  - User must install separately
  - Enable with `-DUSE_HTSLIB` and `-lhts`

**Build-time only**:
- Standard C++ library (libstdc++ or libc++)
- C standard library (glibc or musl)

### Testing Strategy

STAR uses **regression testing** rather than unit tests:

1. **Test Datasets**:
   - Small reference genome (e.g., 1 MB bacterial genome)
   - Synthetic FASTQ reads with known alignments
   - Stored in `test/` or downloaded on demand

2. **Test Commands**:
   ```bash
   # Generate genome index
   STAR --runMode genomeGenerate \
        --genomeDir test/genome \
        --genomeFastaFiles test/chr1.fa

   # Align reads
   STAR --genomeDir test/genome \
        --readFilesIn test/reads.fq \
        --outFileNamePrefix test/output_
   ```

3. **Validation**:
   - Compare SAM output to known-good reference
   - Check mapping statistics in `Log.final.out`
   - Verify splice junction discovery

4. **Continuous Integration**:
   - Not heavily used (as of recent versions)
   - Some projects use TravisCI for basic builds

### Debug Workflow

When developing new features:

1. **Build debug binary**:
   ```bash
   make clean
   make STAR CXXFLAGS="-g -O0 -DDEBUG"
   ```

2. **Run under debugger**:
   ```bash
   gdb --args ./STAR [params]
   (gdb) break alignReads.cpp:123
   (gdb) run
   ```

3. **Enable verbose logging**:
   - Many modules have `#ifdef DEBUG` blocks
   - Add custom logging to trace execution

4. **Test on small dataset**:
   - Use `--genomeLoad NoSharedMemory` (simpler)
   - Single thread (`--runThreadN 1`) avoids race conditions

### Release Process

1. **Version tagging**:
   - Tags in git: `2.7.11a`, `2.7.11b`
   - Version hardcoded in `source/STAR.cpp`

2. **Build matrix**:
   - Linux x86-64 (static binary)
   - macOS Intel (static binary)
   - Source tarball (for custom builds)

3. **Testing**:
   - Run on multiple test datasets
   - Check backward compatibility with old indices
   - Verify against SAMtools validation

4. **Distribution**:
   - GitHub Releases (binary + source)
   - Bioconda package (`conda install star`)
   - Docker images (`quay.io/biocontainers/star`)

### Development Tools

**Recommended**:
- Editor: Any (vim, emacs, VSCode)
- Debugger: GDB (Linux), LLDB (macOS)
- Profiler: `perf`, `gprof`, `Instruments`
- Memory checker: Valgrind, AddressSanitizer

**Code style**:
- No official style guide
- Mix of camelCase and snake_case
- Indentation: 4 spaces
- Comments: C++ style `//`

### Build Performance

**Clean build time** (on modern workstation):
- Full build: ~2-3 minutes (150 source files)
- Incremental rebuild: ~10-30 seconds (touch 1-5 files)
- Parallel make: `make -j8` cuts time in half

**Binary size**:
- Debug build: ~50 MB (with symbols)
- Release build: ~10 MB (unstripped)
- Stripped release: ~2-3 MB

### Common Build Issues

1. **OpenMP not found**:
   - Install `libomp-dev` (Linux) or `libomp` (macOS)
   - Or disable OpenMP (slower): remove `-fopenmp`

2. **zlib missing**:
   - Install `zlib1g-dev` (Debian) or `zlib-devel` (RedHat)

3. **C++11 not supported**:
   - Upgrade compiler: GCC 4.8+, Clang 3.4+

4. **Static linking fails on macOS**:
   - Don't use static build on macOS
   - Use `STARforMacStatic` target instead

## Implications for Modification

When modifying STAR's build and development workflow:

1. **Adding Source Files**: Update Makefile if needed (wildcard should catch new .cpp)
2. **New Dependencies**: Document in README, update Makefile LDFLAGS
3. **Platform Support**: Test on Linux and macOS, verify static builds
4. **Optimization Flags**: Benchmark before/after, test on different CPUs
5. **Testing**: Add regression tests for new features, not just unit tests
6. **Documentation**: Update `doc/STARmanual.pdf` with new parameters
