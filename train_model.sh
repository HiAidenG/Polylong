#bash script for training the model, meant for use on a cluster but you could also adapt this for colab
#WARNING: you must set up your environment beforehand, see README

# scripts from DNABERT repo, link here: https://github.com/jerryji1993/DNABERT

set -euxo pipefail

batch_size = $training_batch

# get the pre-trained model from drive
cd utils 
gdown https://drive.google.com/u/0/uc?id=1BJjqb5Dl2lNMg2warsFQ0-Xvn1xxfFXC&export=download
unzip 6-new-12w-0.zip

cd .. && mkdir model_output && cd DNABERT/examples

python run_finetune.py \
    --model_type dna \
    --tokenizer_name=dna6 \
    --model_name_or_path ../../utils/6-new-12w-0 \
    --task_name dnaprom \
    --do_train \
    --do_eval \
    --data_dir ../../data/507ds/ \
    --max_seq_length 512 \
    --per_gpu_eval_batch_size=$batch_size  \
    --per_gpu_train_batch_size=$batch_size   \
    --learning_rate 2e-4 \
    --num_train_epochs 5.0 \
    --output_dir ../../model_output \
    --evaluate_during_training \
    --logging_steps 100 \
    --save_steps 4000 \
    --warmup_percent 0.1 \
    --hidden_dropout_prob 0.1 \
    --overwrite_output \
    --weight_decay 0.01 \
    --n_process 8

#get predictions on the test set
#there's an issue here with the DNABERT scripts, I'm not sure why you would want predictions on the dev set
#so rather than mess around with their code I just rename the files

mv ../../data/507ds/dev.tsv ../../data/507ds/actual_dev.tsv
mv ../../data/507ds/test.tsv ../../data/507ds/dev.tsv

python run_finetune.py \
    --model_type dna \
    --tokenizer_name=dna6 \
    --model_name_or_path ../../utils/6-new-12w-0 \
    --task_name dnaprom \
    --do_predict \
    --data_dir ../../data/507ds/  \
    --max_seq_length 512 \
    --per_gpu_pred_batch_size=64   \
    --output_dir ../../model_output \
    --predict_dir ../../model_output \
    --n_process 48


