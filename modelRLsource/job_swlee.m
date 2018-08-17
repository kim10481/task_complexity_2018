clear all;
close all;

job_opt.name=14;

list_subject=[0, 1:3 5:12 14:20];
%list_subject=[0];

results=cell(1,length(list_subject));
for i=1:1:length(list_subject)
    sched = findResource('scheduler', 'configuration', 'NeuroEcon.local')
    set(sched,'SubmitArguments', '-l walltime=24:00:00')


    pjob = createParallelJob(sched);
    job_opt.list=list_subject(i);
    set(pjob, 'FileDependencies', {'SIMUL_Arbitration_v3_GA2.m'})
    set(pjob, 'MaximumNumberOfWorkers', 1)
    set(pjob, 'MinimumNumberOfWorkers', 1)
    t = createTask(pjob, @SIMUL_Arbitration_v3_GA2, 1, {job_opt})
    submit(pjob);
    %results{1,i} = getAllOutputArguments(pjob);
end

%errmsgs = get(pjob.Tasks, {'ErrorMessage'});
%nonempty = ~cellfun(@isempty, errmsgs);
%celldisp(errmsgs(nonempty));

%waitForState(pjob);
%taskoutput = getAllOutputArguments(pjob);

disp('### Job submitted. Use "job_msg(id)" for more details  or "job_del(id)" for deleting log files. ###');
