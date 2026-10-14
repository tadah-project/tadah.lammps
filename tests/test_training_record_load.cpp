// pair_style tadah (tadah::LammpsTadah, the object pair_coeff builds from
// the potential file) and the training record of a potential: the comment
// block that tadah train, tadah hpo and tadah refit write above the key
// lines, naming the run's DBFILE, VALIDATION, HPOTARGET, UNCERTAINTY,
// VERBOSE, NUMERIC, OUTFILE, AUDIT, AUDIT_SAMPLE and HEALTH_LOG for the
// user's record (MK's decision of 2026-10-14). The block is comments only:
// pair_coeff reads neither its keys nor its values.

#include "catch2/catch.hpp"

#include <tadah/core/context.h>
#include <tadah/core/logger.h>
#include <tadah/core/periodic_table.h>
#include <tadah/tadah.lammps/lammps_tadah.h>

#include <cstdio>
#include <fstream>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

namespace {

/* A file unique to this process, removed with the object. */
struct TempFile {
  std::string path;
  TempFile(const std::string &stem, const std::string &body)
      : path(stem + "_" + std::to_string(::getpid())) {
    std::ofstream(path) << body;
  }
  ~TempFile() { std::remove(path.c_str()); }
};

/* A Ta pair potential (D2_BP, two functions) as tadah train writes it. */
std::string pair_potential() {
  tadah::core::Context c;
  c.add("ATOM", "Ta");
  c.add("WATOM", 73);
  c.add("BIAS", "false");
  c.add("NORM", "false");
  c.add("INIT2B", "true");
  c.add("INITMB", "false");
  c.add("TYPE2B", "D2_BP 2 2 Ta Ta");
  c.add("CGRID2B", "2.6 3.3");
  c.add("SGRID2B", "1.0 1.0");
  c.add("RCTYPE2B", "Cut_Cos");
  c.add("RCUT2B", 5.0);
  c.add("MODEL", "M_BLR BF_Linear");
  c.add("WEIGHTS", std::vector<double>{0.25, -0.5});
  std::ostringstream os;
  c.write_potential(os, 18);
  return os.str();
}

/* The training record of a run whose configuration set every key a
 * potential keeps only there, VERBOSE 2 among them (a key line VERBOSE 2
 * would set the Tadah log level of the process), with a DBFILE that does
 * not exist. */
const char *const kRecord =
    "# --- training record (not read by tadah) ---\n"
    "# DBFILE no_such_dataset.db other.db\n"
    "# VALIDATION config.val\n"
    "# HPOTARGET config.hpo\n"
    "# UNCERTAINTY true\n"
    "# VERBOSE 2\n"
    "# NUMERIC 9\n"
    "# OUTFILE pot.tadah\n"
    "# AUDIT error\n"
    "# AUDIT_SAMPLE 64\n"
    "# HEALTH_LOG full\n"
    "# --- end of training record ---\n";

/* pair_coeff * * <file> Ta: LammpsTadah from the file, then the type map. */
std::unique_ptr<tadah::LammpsTadah> load(const std::string &path) {
  std::vector<std::string> args = {"*", "*", path, "Ta"};
  std::vector<char *> argv;
  for (auto &a : args) argv.push_back(&a[0]);
  std::unique_ptr<tadah::LammpsTadah> lt(
      new tadah::LammpsTadah(static_cast<int>(argv.size()), argv.data()));
  lt->initialize(static_cast<int>(argv.size()), argv.data());
  return lt;
}

}  // namespace

/* Pins: pair_coeff loads a potential whose key lines follow its training
 * record (DBFILE of a dataset that does not exist, VALIDATION, HPOTARGET,
 * UNCERTAINTY true, VERBOSE 2, NUMERIC 9, OUTFILE, AUDIT error,
 * AUDIT_SAMPLE 64, HEALTH_LOG full) to the same context, key by key, as the
 * same key lines without the block (a potential without a record, as every
 * earlier version wrote, loads as before), with none of the ten keys but
 * VERBOSE, which the model's construction stores at its default as for the
 * potential without the block (not the record's 2), and leaves the Tadah
 * log level where the potential without the block leaves it.
 * Before the fix: passes as well (a pin: core::CoreUserConfigParser, the
 * reader of pair_coeff, skips every line whose first non-blank character
 * is '#', and the training record is a block of such lines). */
TEST_CASE("pair_style tadah ignores the training record of a potential",
          "[lammps_tadah][training_record]") {
  tadah::core::PeriodicTable::initialize();
  auto &log = tadah::core::Logger::instance();
  const auto old_level = log.getLogLevel();
  const std::string body = pair_potential();
  const TempFile bare("training_record_bare.tadah", body);
  const TempFile with("training_record_with.tadah", kRecord + body);

  log.setLogLevel(tadah::core::LogLevel::WARNING);
  std::unique_ptr<tadah::LammpsTadah> lt_bare = load(bare.path);
  const auto level_bare = log.getLogLevel();
  log.setLogLevel(tadah::core::LogLevel::WARNING);
  std::unique_ptr<tadah::LammpsTadah> lt_with = load(with.path);
  const auto level_with = log.getLogLevel();
  log.setLogLevel(old_level);

  REQUIRE(lt_bare != nullptr);
  REQUIRE(lt_with != nullptr);
  CHECK(lt_with->c.storage_snapshot() == lt_bare->c.storage_snapshot());
  for (const char *key : {"DBFILE", "VALIDATION", "HPOTARGET", "UNCERTAINTY",
                          "NUMERIC", "OUTFILE", "AUDIT", "AUDIT_SAMPLE",
                          "HEALTH_LOG"}) {
    INFO("key " << key);
    CHECK_FALSE(lt_with->c.has_stored(key));
  }
  // The model's construction stores VERBOSE at its blueprint default, as
  // for the potential without the block: not the record's 2.
  CHECK(lt_with->c.stored_value("VERBOSE") ==
        lt_bare->c.stored_value("VERBOSE"));
  CHECK(lt_with->c.stored_value("VERBOSE") != std::optional<std::string>("2"));
  CHECK(level_with == level_bare);
}
