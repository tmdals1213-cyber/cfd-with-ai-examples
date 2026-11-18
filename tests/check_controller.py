import importlib.util,tempfile,pathlib,unittest
ROOT=pathlib.Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('controller',ROOT/'scripts/case.py');c=importlib.util.module_from_spec(spec);spec.loader.exec_module(c)
class AuditTests(unittest.TestCase):
 def fixture(self):
  t=tempfile.TemporaryDirectory();self.addCleanup(t.cleanup);p=pathlib.Path(t.name);(p/'run.log').write_text('normal solver end');return p
 def test_os_zero_abort_is_failure(self):
  p=self.fixture();(p/'run.log').write_text('MPI_ABORT was invoked with code zero')
  with self.assertRaises(AssertionError):c.audit(p,.02)
 def test_finite_finished_custom_history(self):
  p=self.fixture();(p/'history.csv').write_text('step,time_s,total_energy_error_J,div_minus_source_s_1,pressure_relative_error,mass_balance_relative,phase_power_W\n1,.02,1e-9,1e-8,8e-9,1e-8,1\n')
  self.assertEqual(c.audit(p,.02)['rows'],1)
 def test_incomplete_time_refuses(self):
  p=self.fixture();(p/'history.csv').write_text('step,time_s\n1,.01\n')
  with self.assertRaises(AssertionError):c.audit(p,.02)
 def test_excess_mass_error_refuses(self):
  p=self.fixture();(p/'history.csv').write_text('step,time_s,mass_balance_relative,phase_power_W\n1,.02,.01,1\n')
  with self.assertRaises(AssertionError):c.audit(p,.02)
if __name__=='__main__':unittest.main()
